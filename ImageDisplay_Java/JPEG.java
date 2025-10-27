import java.awt.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;
import java.util.ArrayList;
import java.util.List;
import static java.lang.System.exit;

public class JPEG {

    JFrame frame;
    JLabel lbIm1;
    JLabel lbIm2; // Label for the output image
    BufferedImage imgOne;
    BufferedImage imgTwo; // BufferedImage for the output image

    int width = 352;
    int height = 288;

    // Placeholder for storing quantized DCT coefficients for all blocks and channels
    double[][][][][] allCoefficients = new double[3][height / 8][width / 8][8][8]; // [Channel][BlockY][BlockX][CoeffIndex]

    /**
     * Represents a chunk of data for GUI update (e.g., a decoded block).
     * We'll send the block coordinates and the 8x8 pixel data.
     */
    private static class DecodeUpdate {
        int blockX, blockY, channel;
        int[][] pixels = new int[8][8]; // Store final RGB pixel values for the block

        DecodeUpdate(int ch, int by, int bx, int[][] blockPixels) {
            this.channel = ch;
            this.blockY = by;
            this.blockX = bx;
            // Copy the pixel data
            for(int y=0; y<8; y++){
                System.arraycopy(blockPixels[y], 0, this.pixels[y], 0, 8);
            }
        }
    }

    /**
     * SwingWorker for handling the decoding simulation in the background.
     */
    private class DecoderWorker extends SwingWorker<Void, DecodeUpdate> {
        private final int quantizationLevelN;
        private final int deliveryModeM;
        private final int latencyL;

        DecoderWorker(int n, int m, int l) {
            this.quantizationLevelN = n;
            this.deliveryModeM = m;
            this.latencyL = l;
        }

        private double[][] bitmask(int bitShift, double[][] original){
            double[][] ret = new double[original.length][original[0].length];
            for (int i = 0; i < original.length; i++) {
                for (int i1 = 0; i1 < original[i].length; i1++) {
                    double original_val = original[i][i1];
                    int magnitude = (int)Math.abs(Math.round(original_val));
                    magnitude>>=bitShift;
                    magnitude<<=bitShift;
                    System.out.println("magnitude is " + magnitude);
                    if(original_val < 0){
                        ret[i][i1] = -magnitude;
                    }else{
                        ret[i][i1] = magnitude;
                    }
                }
            }
            return ret;
        }

        @Override
        protected Void doInBackground() throws Exception {
            // Initialize imgTwo (e.g., make it black or grey initially)
            Graphics2D g2d = imgTwo.createGraphics();
            g2d.setColor(Color.LIGHT_GRAY); // Or Color.BLACK
            g2d.fillRect(0, 0, width, height);
            g2d.dispose();
            publish(); // Publish initial blank image state (optional)

            if (deliveryModeM == 1) {
                // Mode 1: Baseline (Sequential)
                for (int blockY = 0; blockY < height / 8; blockY++) {
                    for (int blockX = 0; blockX < width / 8; blockX++) {
                        int[][] blockPixelsRGB = new int[8][8]; // Store final RGB for this block

                        // Process all 3 channels for this block before updating GUI
                        for (int channel = 0; channel < 3; channel++) {
                            // 5. Dequantize
                            double[][] quantizedCoeffs = allCoefficients[channel][blockY][blockX];
                            double[][] dequantizedCoeffs = dequantize(quantizedCoeffs, quantizationLevelN);

                            // 6. Inverse DCT
                            double[][] reconstructedBlock = inverseDCT(dequantizedCoeffs);

                            // Combine channels into blockPixelsRGB
                            for (int y = 0; y < 8; y++) {
                                for (int x = 0; x < 8; x++) {
                                    double shiftedValue = reconstructedBlock[y][x] + 128.0;
                                    int pixelCompValue = (int) Math.round(Math.max(0.0, Math.min(255.0, shiftedValue)));

                                    int existingRgb = (channel == 0) ? 0xFF000000 : blockPixelsRGB[y][x]; // Start fresh for R, reuse for G, B
                                    int r = (existingRgb >> 16) & 0xFF;
                                    int g = (existingRgb >> 8) & 0xFF;
                                    int b = existingRgb & 0xFF;

                                    if (channel == 0) r = pixelCompValue;
                                    else if (channel == 1) g = pixelCompValue;
                                    else b = pixelCompValue;

                                    blockPixelsRGB[y][x] = 0xFF000000 | (r << 16) | (g << 8) | b;
                                }
                            }
                        } // End channel loop for this block

                        // Publish a single RGB block
                        publish(new DecodeUpdate(-1, blockY, blockX, blockPixelsRGB)); // Use channel=-1 to signify RGB block

                        // Sleep after processing one block
                        if (latencyL > 0) {
                            Thread.sleep(latencyL);
                        }
                    } // End blockX loop
                } // End blockY loop

            }
            else if (deliveryModeM == 2) {
                //zigzag traversal
                int i, j;
                double[][][][][] spectrumCoefficients = new double[3][height / 8][width / 8][8][8];
                int total = 0;
                while(total <= 14){
                    List<DecodeUpdate> updates = new ArrayList<>();
                    for(i = 0, j = total - i; j>=0; i++, j--){
                        if(i > 7 || j > 7) continue;
                        for (int blockY = 0; blockY < height / 8; blockY++) {
                            for (int blockX = 0; blockX < width / 8; blockX++) {
                                int[][] blockPixelsRGB = new int[8][8]; // Store final RGB for this block

                                // Process all 3 channels for this block before updating GUI
                                for (int channel = 0; channel < 3; channel++) {
                                    spectrumCoefficients[channel][blockY][blockX][i][j] = allCoefficients[channel][blockY][blockX][i][j];
                                    // 5. Dequantize
                                    double[][] quantizedCoeffs = spectrumCoefficients[channel][blockY][blockX];
                                    double[][] dequantizedCoeffs = dequantize(quantizedCoeffs, quantizationLevelN);

                                    // 6. Inverse DCT
                                    double[][] reconstructedBlock = inverseDCT(dequantizedCoeffs);

                                    // Combine channels into blockPixelsRGB
                                    for (int y = 0; y < 8; y++) {
                                        for (int x = 0; x < 8; x++) {
                                            double shiftedValue = reconstructedBlock[y][x] + 128.0;
                                            int pixelCompValue = (int) Math.round(Math.max(0.0, Math.min(255.0, shiftedValue)));

                                            int existingRgb = (channel == 0) ? 0xFF000000 : blockPixelsRGB[y][x]; // Start fresh for R, reuse for G, B
                                            int r = (existingRgb >> 16) & 0xFF;
                                            int g = (existingRgb >> 8) & 0xFF;
                                            int b = existingRgb & 0xFF;

                                            if (channel == 0) r = pixelCompValue;
                                            else if (channel == 1) g = pixelCompValue;
                                            else b = pixelCompValue;

                                            blockPixelsRGB[y][x] = 0xFF000000 | (r << 16) | (g << 8) | b;
                                        }
                                    }
                                } // End channel loop for this block
                                updates.add(new DecodeUpdate(-1, blockY, blockX, blockPixelsRGB));

                            } // End blockX loop
                        } // End blockY loop
                    }
                    total++;
                    // Publish all blocks with selected spectrum
                    publish(updates.toArray(new DecodeUpdate[0]));
                    // Sleep after processing one block
                    if (latencyL > 0) {
                        Thread.sleep(latencyL);
                    }
                }

            }
            else if (deliveryModeM == 3) {
                System.out.println("in mode 3");
                int shifted = 7;
                for (int i = shifted; i >=1; i--) {
                    List<DecodeUpdate> updates = new ArrayList<>();
                    System.out.println("shifted is " + shifted);
                    for (int blockY = 0; blockY < height / 8; blockY++) {
                        for (int blockX = 0; blockX < width / 8; blockX++) {
                            int[][] blockPixelsRGB = new int[8][8]; // Store final RGB for this block

                            // Process all 3 channels for this block before updating GUI
                            for (int channel = 0; channel < 3; channel++) {
                                // 5. Dequantize
                                double[][] maskedQuantizedCoeffs = bitmask(i, allCoefficients[channel][blockY][blockX]);

                                double[][] dequantizedCoeffs = dequantize(maskedQuantizedCoeffs, quantizationLevelN);

                                // 6. Inverse DCT
                                double[][] reconstructedBlock = inverseDCT(dequantizedCoeffs);
                                // Combine channels into blockPixelsRGB
                                for (int y = 0; y < 8; y++) {
                                    for (int x = 0; x < 8; x++) {
                                        double shiftedValue = reconstructedBlock[y][x] + 128.0;
                                        int pixelCompValue = (int) Math.round(Math.max(0.0, Math.min(255.0, shiftedValue)));

                                        int existingRgb = (channel == 0) ? 0xFF000000 : blockPixelsRGB[y][x]; // Start fresh for R, reuse for G, B
                                        int r = (existingRgb >> 16) & 0xFF;
                                        int g = (existingRgb >> 8) & 0xFF;
                                        int b = existingRgb & 0xFF;

                                        if (channel == 0) r = pixelCompValue;
                                        else if (channel == 1) g = pixelCompValue;
                                        else b = pixelCompValue;

                                        blockPixelsRGB[y][x] = 0xFF000000 | (r << 16) | (g << 8) | b;
                                    }
                                }
                            } // End channel loop for this block

                            updates.add(new DecodeUpdate(-1, blockY, blockX, blockPixelsRGB));
                        } // End blockX loop
                    } // End blockY loop
                    // Publish all with bit masked
                    publish(updates.toArray(new DecodeUpdate[0]));
                    // Sleep after processing one block
                    if (latencyL > 0) {
                        Thread.sleep(latencyL);
                    }
                }

            }

            System.out.println("DecoderWorker finished background task.");
            return null;
        }
        @Override
        protected void process(List<DecodeUpdate> chunks) {
            // This runs on the GUI thread
            for (DecodeUpdate update : chunks) {// in baseline mode, one block one be rendered one time; in other modes, the whole pic will be rendered, though diff resolution
                    for (int y = 0; y < 8; y++) {
                        for (int x = 0; x < 8; x++) {
                            int absoluteY = update.blockY * 8 + y;
                            int absoluteX = update.blockX * 8 + x;
                            if (absoluteX < width && absoluteY < height) { // Bounds check
                                imgTwo.setRGB(absoluteX, absoluteY, update.pixels[y][x]);
                            }
                        }
                    }
            }
            if (lbIm2 != null) {
                lbIm2.setIcon(new ImageIcon(imgTwo)); // Update the label's image icon
                // frame.pack(); // Might cause resizing, repaint is usually enough
                frame.repaint();
            }
        }

        @Override
        protected void done() {
            // Called on the GUI thread when doInBackground finishes
            System.out.println("DecoderWorker is done.");
            // You could re-enable buttons or show a completion message here
        }
    }

    /**
     * Reads the raw RGB image file.
     * Assumes the format is RRR...GGG...BBB...
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
            System.err.println("Error reading image file: " + e.getMessage());
            exit(-1);
        }
    }

    /**
     * The main processing function.
     * Orchestrates encoding and decoding based on input parameters.
     */
    public void process(String[] args){
        String imagePath = args[0];
        int quantizationLevelN = Integer.parseInt(args[1]); // N
        int deliveryModeM = Integer.parseInt(args[2]);      // M
        int latencyL = Integer.parseInt(args[3]);           // L

        // 1. Read Input Image
        imgOne = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        readImageRGB(width, height, imagePath, imgOne);

        // Initialize output image (can be initially black or a copy)
        imgTwo = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        // --- Start of ENCODER Logic ---
        for (int channel = 0; channel < 3; channel++) { // 0=R, 1=G, 2=B
            for (int blockY = 0; blockY < height / 8; blockY++) {
                for (int blockX = 0; blockX < width / 8; blockX++) {

                    // 2. Break channel into 8x8 blocks
                    double[][] currentBlock = new double[8][8];
                    for (int y = 0; y < 8; y++) {
                        for (int x = 0; x < 8; x++) {
                            int absoluteY = blockY * 8 + y;
                            int absoluteX = blockX * 8 + x;
                            int rgb = imgOne.getRGB(absoluteX, absoluteY);
                            int componentValue;
                            if (channel == 0) componentValue = (rgb >> 16) & 0xFF; // R
                            else if (channel == 1) componentValue = (rgb >> 8) & 0xFF;  // G
                            else componentValue = rgb & 0xFF;                     // B
                            currentBlock[y][x] = componentValue - 128.0;
                        }
                    }

                    // 3. Discrete Cosine Transform
                    double[][] dctCoefficients = forwardDCT(currentBlock);
                    // 4. Quantize based on table (uniform 2^N)
                    allCoefficients[channel][blockY][blockX] = quantize(dctCoefficients, quantizationLevelN);
                }
            }
        }
        // --- Start of DECODER Logic (Simulation) ---
        DecoderWorker decoderWorker = new DecoderWorker(quantizationLevelN, deliveryModeM, latencyL);
        decoderWorker.execute();
        displayImages(true); // Temporarily show GUI instantly
        // --- End of DECODER Logic ---
    }


    /**
     * Sets up and displays the GUI frame with the two images.
     * @param makeVisible If true, makes the frame visible.
     */
    private void displayImages(boolean makeVisible) {
        if (frame == null) {
            frame = new JFrame();
            GridBagLayout gLayout = new GridBagLayout();
            frame.getContentPane().setLayout(gLayout);

            // Use imgOne for both initially, imgTwo will be updated by decoder
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
            c.gridy = 0;
            frame.getContentPane().add(lbIm2, c);

            frame.pack();
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        }
        if (makeVisible) {
            frame.setVisible(true);
        }
    }

    // Helper method to create a deep copy of a BufferedImage
    private BufferedImage deepCopy(BufferedImage bi) {
        ColorModel cm = bi.getColorModel();
        boolean isAlphaPremultiplied = cm.isAlphaPremultiplied();
        WritableRaster raster = bi.copyData(bi.getRaster().createCompatibleWritableRaster());
        return new BufferedImage(cm, raster, isAlphaPremultiplied, null);
    }

    /**
     * Calculates the normalization coefficient C(k) for DCT/IDCT.
     * @param k The frequency index (u or v).
     * @return The normalization coefficient.
     */
    private double C(int k) {
        if (k == 0) {
            return 1.0 / Math.sqrt(2.0);
        } else {
            return 1.0;
        }
    }

    /**
     * Performs the 2D Forward Discrete Cosine Transform on an 8x8 block.
     * Input block values should be in the range [-128, 127].
     * @param block The 8x8 input block f(x,y).
     * @return The 8x8 block of DCT coefficients F(u,v).
     */
    public double[][] forwardDCT(double[][] block) {
        double[][] coefficients = new double[8][8];
        double pi = Math.PI;

        for (int u = 0; u < 8; u++) {
            for (int v = 0; v < 8; v++) {
                double sum = 0.0;
                for (int x = 0; x < 8; x++) {
                    for (int y = 0; y < 8; y++) {
                        double cosX = Math.cos(((2.0 * x + 1.0) * u * pi) / 16.0);
                        double cosY = Math.cos(((2.0 * y + 1.0) * v * pi) / 16.0);
                        sum += block[y][x] * cosX * cosY;
                    }
                }
                coefficients[v][u] = (1.0 / 4.0) * C(u) * C(v) * sum;
            }
        }
        return coefficients;
    }

    /**
     * Performs the 2D Inverse Discrete Cosine Transform on an 8x8 block of coefficients.
     * Output block values will be approximately in the range [-128, 127].
     * @param coefficients The 8x8 block of DCT coefficients F(u,v).
     * @return The reconstructed 8x8 block f(x,y).
     */
    public double[][] inverseDCT(double[][] coefficients) {
        double[][] block = new double[8][8];
        double pi = Math.PI;

        for (int x = 0; x < 8; x++) {
            for (int y = 0; y < 8; y++) {
                double sum = 0.0;
                for (int u = 0; u < 8; u++) {
                    for (int v = 0; v < 8; v++) {
                        double cosX = Math.cos(((2.0 * x + 1.0) * u * pi) / 16.0);
                        double cosY = Math.cos(((2.0 * y + 1.0) * v * pi) / 16.0);
                        sum += C(u) * C(v) * coefficients[v][u] * cosX * cosY;
                    }
                }
                block[y][x] = (1.0 / 4.0) * sum;
            }
        }
        return block;
    }

    /**
     * Quantizes the DCT coefficients using a uniform step size of 2^N.
     * @param coeffs The 8x8 block of DCT coefficients F(u,v).
     * @param n The quantization level (0-7).
     * @return The 8x8 block of quantized coefficients F'(u,v).
     */
    public double[][] quantize(double[][] coeffs, int n) {
        if (n < 0 || n > 7) {
            System.err.println("Warning: Quantization level N must be between 0 and 7. Clamping.");
            n = Math.max(0, Math.min(7, n));
        }

        double[][] quantizedCoeffs = new double[8][8];
        double quantizationStep = Math.pow(2, n);

        // Handle N=0 case separately to avoid division by 1 unnecessarily
        if (n == 0) {
            // No quantization, just copy (or return original)
            for (int v = 0; v < 8; v++) {
                System.arraycopy(coeffs[v], 0, quantizedCoeffs[v], 0, 8);
            }
            return quantizedCoeffs;
        }

        for (int v = 0; v < 8; v++) { // Row index (maps to v in F(u,v))
            for (int u = 0; u < 8; u++) { // Column index (maps to u in F(u,v))
                quantizedCoeffs[v][u] = Math.round(coeffs[v][u] / quantizationStep);
            }
        }
        return quantizedCoeffs;
    }

    /**
     * Dequantizes the DCT coefficients using a uniform step size of 2^N.
     * @param quantizedCoeffs The 8x8 block of quantized coefficients F'(u,v).
     * @param n The quantization level (0-7).
     * @return The 8x8 block of dequantized coefficients F(u,v) (approximate).
     */
    public double[][] dequantize(double[][] quantizedCoeffs, int n) {
        if (n < 0 || n > 7) {
            System.err.println("Warning: Quantization level N must be between 0 and 7. Clamping.");
            n = Math.max(0, Math.min(7, n));
        }

        double[][] coeffs = new double[8][8];
        double quantizationStep = Math.pow(2, n);

        // Handle N=0 case separately
        if (n == 0) {
            for (int v = 0; v < 8; v++) {
                System.arraycopy(quantizedCoeffs[v], 0, coeffs[v], 0, 8);
            }
            return coeffs;
        }

        for (int v = 0; v < 8; v++) {
            for (int u = 0; u < 8; u++) {
                coeffs[v][u] = quantizedCoeffs[v][u] * quantizationStep;
            }
        }
        return coeffs;
    }


    public static void main(String[] args) {
        if (args.length < 4) {
            System.err.println("Usage: java JPEG <InputImage> <QuantizationLevel N> <DeliveryMode M> <Latency L>");
            System.err.println("N: 0-7");
            System.err.println("M: 1 (Baseline), 2 (Spectral), 3 (Successive Bit)");
            System.err.println("L: Latency in milliseconds (e.g., 0, 100)");
            return;
        }

        // Basic argument validation
        try {
            Integer.parseInt(args[1]); // N
            int mode = Integer.parseInt(args[2]); // M
            if (mode < 1 || mode > 3) throw new NumberFormatException();
            Integer.parseInt(args[3]); // L
        } catch (NumberFormatException e) {
            System.err.println("Invalid number format for N, M, or L.");
            System.err.println("Usage: java JPEG <InputImage> <QuantizationLevel N> <DeliveryMode M> <Latency L>");
            return;
        }

        JPEG jpegProcessor = new JPEG();
        // Use SwingUtilities.invokeLater to ensure GUI setup happens on the Event Dispatch Thread
        SwingUtilities.invokeLater(() -> jpegProcessor.process(args));
    }
}