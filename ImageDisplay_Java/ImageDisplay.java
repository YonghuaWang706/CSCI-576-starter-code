import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;

public class ImageDisplay {

    JFrame frame;
    JLabel lbIm1;
    JLabel lbIm2; // Label for the output image
    BufferedImage imgOne;
    BufferedImage imgTwo; // BufferedImage for the output image

//    int width = 512;
//    int height = 512;
    int width = 352;
    int height = 288;

    // Custom inner class to hold YUV data, using doubles for precision
    // and to allow for negative U and V values.
    public static class YUVImage {
        double[][] y;
        double[][] u;
        double[][] v;

        public YUVImage(int width, int height) {
            y = new double[height][width];
            u = new double[height][width];
            v = new double[height][width];
        }
    }

    /**
     * Calculates the sum of absolute differences between two RGB images.
     * @param original The original source image.
     * @param quantized The processed, quantized image.
     * @return The total error as a long integer.
     */
    private long calculateError(BufferedImage original, BufferedImage quantized) {
        long totalError = 0;
        int width = original.getWidth();
        int height = original.getHeight();

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int rgbOrig = original.getRGB(x, y);
                int rOrig = (rgbOrig >> 16) & 0xFF;
                int gOrig = (rgbOrig >> 8) & 0xFF;
                int bOrig = rgbOrig & 0xFF;

                int rgbQuant = quantized.getRGB(x, y);
                int rQuant = (rgbQuant >> 16) & 0xFF;
                int gQuant = (rgbQuant >> 8) & 0xFF;
                int bQuant = rgbQuant & 0xFF;

                totalError += Math.abs(rOrig - rQuant);
                totalError += Math.abs(gOrig - gQuant);
                totalError += Math.abs(bOrig - bQuant);
            }
        }
        return totalError;
    }

    /**
     * Converts a BufferedImage from RGB color space to a YUVImage object.
     */
    private YUVImage convertToYUV(BufferedImage img) {
        int width = img.getWidth();
        int height = img.getHeight();
        YUVImage yuvImage = new YUVImage(width, height);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int rgb = img.getRGB(x, y);
                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = rgb & 0xFF;

                yuvImage.y[y][x] = 0.299 * r + 0.587 * g + 0.114 * b;
                yuvImage.u[y][x] = -0.14713 * r - 0.28886 * g + 0.436 * b;
                yuvImage.v[y][x] = 0.615 * r - 0.51499 * g - 0.10001 * b;
            }
        }
        return yuvImage;
    }

    /**
     * Converts a YUVImage object back to a BufferedImage in RGB format.
     */
    private BufferedImage convertToRGB(YUVImage yuvImage) {
        int width = yuvImage.y[0].length;
        int height = yuvImage.y.length;
        BufferedImage rgbImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double yVal = yuvImage.y[y][x];
                double uVal = yuvImage.u[y][x];
                double vVal = yuvImage.v[y][x];

                double r_double = yVal + 1.13983 * vVal;
                double g_double = yVal - 0.39465 * uVal - 0.58060 * vVal;
                double b_double = yVal + 2.03211 * uVal;

                int r = (int) Math.max(0, Math.min(255, r_double));
                int g = (int) Math.max(0, Math.min(255, g_double));
                int b = (int) Math.max(0, Math.min(255, b_double));

                int pix = 0xFF000000 | (r << 16) | (g << 8) | b;
                rgbImage.setRGB(x, y, pix);
            }
        }
        return rgbImage;
    }

    /**
     * Performs uniform quantization on a single color channel.
     */
    private void uniformQuantizeChannel(double[][] channel, int numBits, double minRange, double maxRange) {
        if (numBits > 8) return; // Ignore invalid bit counts

        // BUG FIX: Handle numBits = 0 case
        if (numBits == 0) {
            double neutralValue = (minRange + maxRange) / 2.0;
            if (minRange < 0) neutralValue = 0; // For U/V, neutral is 0
            for (int y = 0; y < channel.length; y++) {
                for (int x = 0; x < channel[0].length; x++) {
                    channel[y][x] = neutralValue;
                }
            }
            return;
        }
        if (numBits < 1) return; // Should not happen with the zero check, but good practice

        int levels = 1 << numBits;
        double step = (maxRange - minRange) / levels;

        for (int y = 0; y < channel.length; y++) {
            for (int x = 0; x < channel[0].length; x++) {
                double originalValue = channel[y][x];
                int levelIndex = (int) Math.floor((originalValue - minRange) / step);
                levelIndex = Math.max(0, Math.min(levels - 1, levelIndex));
                double representativeValue = minRange + (levelIndex + 0.5) * step;
                channel[y][x] = representativeValue;
            }
        }
    }

    /**
     * Performs non-uniform "smart" quantization on a channel using a histogram.
     */
    private void smartQuantizeChannel(double[][] channel, int numBits, int minRange, int maxRange) {
        if (numBits > 8) return;

        // BUG FIX: Handle numBits = 0 case
        if (numBits == 0) {
            double neutralValue = (minRange + maxRange) / 2.0;
            if (minRange < 0) neutralValue = 0; // For U/V, neutral is 0
            for (int y = 0; y < channel.length; y++) {
                for (int x = 0; x < channel[0].length; x++) {
                    channel[y][x] = neutralValue;
                }
            }
            return;
        }
        if (numBits < 1) return;

        int rangeSize = maxRange - minRange + 1;
        int offset = -minRange;
        int[] histogram = new int[rangeSize];
        int totalPixels = channel.length * channel[0].length;

        // 1. Build Histogram
        for (int y = 0; y < channel.length; y++) {
            for (int x = 0; x < channel[0].length; x++) {
                int val = (int) Math.round(channel[y][x]);
                int clampedVal = Math.max(minRange, Math.min(maxRange, val));
                histogram[clampedVal + offset]++;
            }
        }

        int levels = 1 << numBits;
        int pixelsPerBin = totalPixels / levels;
        int[] boundaries = new int[levels + 1];
        boundaries[0] = 0;
        boundaries[levels] = rangeSize - 1;

        // 2. Determine Boundaries based on equal pixel count per bin
        int pixelCount = 0;
        int binIndex = 1;
        for (int i = 0; i < rangeSize; i++) {
            pixelCount += histogram[i];
            if (pixelCount >= pixelsPerBin && binIndex < levels) {
                boundaries[binIndex] = i;
                binIndex++;
                pixelCount = 0; // Reset for next bin
            }
        }

        // 3. Calculate Representative Centroids
        double[] representatives = new double[levels];
        for (int i = 0; i < levels; i++) {
            long valueSum = 0;
            int countInBin = 0;
            for (int j = boundaries[i]; j <= boundaries[i+1]; j++) {
                valueSum += (long) (j - offset) * histogram[j];
                countInBin += histogram[j];
            }
            if (countInBin > 0) {
                representatives[i] = (double) valueSum / countInBin;
            } else {
                // Handle empty bin case - use midpoint of boundary
                representatives[i] = ((boundaries[i] + boundaries[i+1]) / 2.0) - offset;
            }
        }

        // 4. Create Lookup Table and 5. Apply Mapping
        double[] lookupTable = new double[rangeSize];
        for (int i = 0; i < levels; i++) {
            for (int j = boundaries[i]; j <= boundaries[i+1]; j++) {
                lookupTable[j] = representatives[i];
            }
        }

        for (int y = 0; y < channel.length; y++) {
            for (int x = 0; x < channel[0].length; x++) {
                int val = (int) Math.round(channel[y][x]);
                int clampedVal = Math.max(minRange, Math.min(maxRange, val));
                channel[y][x] = lookupTable[clampedVal + offset];
            }
        }
    }

    /**
     * Reads the raw RGB image file.
     */
    private void readImageRGB(int width, int height, String imgPath, BufferedImage img) {
        try {
            int frameLength = width * height * 3;
            File file = new File(imgPath);
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            raf.seek(0);
            byte[] bytes = new byte[frameLength];
            raf.read(bytes);
            raf.close();

            int ind = 0;
            for(int y = 0; y < height; y++) {
                for(int x = 0; x < width; x++) {
                    byte r = bytes[ind];
                    byte g = bytes[ind + height * width];
                    byte b = bytes[ind + height * width * 2];
                    int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                    img.setRGB(x, y, pix);
                    ind++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void showIms(String[] args){
        String imagePath = args[0];
        int colorMode = Integer.parseInt(args[1]);
        int quantizationMode = Integer.parseInt(args[2]);
        int q1 = Integer.parseInt(args[3]);
        int q2 = Integer.parseInt(args[4]);
        int q3 = Integer.parseInt(args[5]);

        imgOne = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        readImageRGB(width, height, imagePath, imgOne);
        imgTwo = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        if (colorMode == 1) { // RGB Space
            YUVImage rgbAsYUV = new YUVImage(width, height);
            for(int y=0; y<height; y++) {
                for(int x=0; x<width; x++) {
                    int rgb = imgOne.getRGB(x,y);
                    rgbAsYUV.y[y][x] = (rgb >> 16) & 0xFF; // R
                    rgbAsYUV.u[y][x] = (rgb >> 8) & 0xFF;  // G
                    rgbAsYUV.v[y][x] = rgb & 0xFF;         // B
                }
            }

            if (quantizationMode == 1) {
                uniformQuantizeChannel(rgbAsYUV.y, q1, 0, 255);
                uniformQuantizeChannel(rgbAsYUV.u, q2, 0, 255);
                uniformQuantizeChannel(rgbAsYUV.v, q3, 0, 255);
            } else { // M=2
                smartQuantizeChannel(rgbAsYUV.y, q1, 0, 255);
                smartQuantizeChannel(rgbAsYUV.u, q2, 0, 255);
                smartQuantizeChannel(rgbAsYUV.v, q3, 0, 255);
            }

            for(int y=0; y<height; y++) {
                for(int x=0; x<width; x++) {
                    int r = (int)Math.round(rgbAsYUV.y[y][x]);
                    int g = (int)Math.round(rgbAsYUV.u[y][x]);
                    int b = (int)Math.round(rgbAsYUV.v[y][x]);
                    int pix = 0xFF000000 | (r << 16) | (g << 8) | b;
                    imgTwo.setRGB(x,y,pix);
                }
            }
        } else if (colorMode == 2) { // YUV Space
            YUVImage yuvImage = convertToYUV(imgOne);
            if (quantizationMode == 1) {
                uniformQuantizeChannel(yuvImage.y, q1, 0, 255);
                uniformQuantizeChannel(yuvImage.u, q2, -112, 112);
                uniformQuantizeChannel(yuvImage.v, q3, -157, 157);
            } else { // M=2
                smartQuantizeChannel(yuvImage.y, q1, 0, 255);
                smartQuantizeChannel(yuvImage.u, q2, -112, 112);
                smartQuantizeChannel(yuvImage.v, q3, -157, 157);
            }
            imgTwo = convertToRGB(yuvImage);
        }

        // Calculate and print the error
        long error = calculateError(imgOne, imgTwo);
        System.out.println("Parameters: <C=" + colorMode + ", M=" + quantizationMode + ", Q1=" + q1 + ", Q2=" + q2 + ", Q3=" + q3 + ">");
        System.out.println("Total Error: " + error);


        frame = new JFrame();
        GridBagLayout gLayout = new GridBagLayout();
        frame.getContentPane().setLayout(gLayout);
        lbIm1 = new JLabel(new ImageIcon(imgOne));
        lbIm2 = new JLabel(new ImageIcon(imgTwo));
        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.weightx = 0.5;
        c.gridx = 0;
        c.gridy = 0;
        frame.getContentPane().add(lbIm1, c);
        c.gridx = 1;
        frame.getContentPane().add(lbIm2, c);
        frame.pack();
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    public static void main(String[] args) {
        if (args.length < 6) {
            System.out.println("Usage: java ImageDisplay <image_path> <C> <M> <Q1> <Q2> <Q3>");
            return;
        }
        ImageDisplay ren = new ImageDisplay();
        ren.showIms(args);
    }
}
