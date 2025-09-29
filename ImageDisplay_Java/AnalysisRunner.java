import java.util.ArrayList;
import java.util.List;
public class AnalysisRunner {
    public static void main(String[] args) {
        if (args.length < 1) {
            System.out.println("Usage: java AnalysisRunner <image_path.rgb>");
            return;
        }

        String imagePath = args[0];
        int[] N_values = {4, 6, 8};

        // Generate commands for each N
        for (int N : N_values) {
            System.out.println("\n// Commands for N = " + N);
            List<int[]> combinations = findCombinations(N);

            for (int[] combo : combinations) {
                int q1 = combo[0];
                int q2 = combo[1];
                int q3 = combo[2];

                // Loop through Color Modes (C=1 for RGB, C=2 for YUV)
                for (int c = 1; c <= 2; c++) {
                    // Loop through Quantization Modes (M=1 for Uniform, M=2 for Smart)
                    for (int m = 1; m <= 2; m++) {
                        // Print the command to be executed
                        System.out.printf("java ImageDisplay %s %d %d %d %d %d --no-gui\n",
                                imagePath, c, m, q1, q2, q3);
                    }
                }
            }
        }
    }

    /**
     * Finds all unique combinations of three non-negative integers (q1, q2, q3)
     * that sum up to a target value N. Each integer can be at most 8.
     * @param N The target sum.
     * @return A list of integer arrays, where each array is a valid combination.
     */
    private static List<int[]> findCombinations(int N) {
        List<int[]> result = new ArrayList<>();
        int maxVal = 8; // Max bits per channel

        for (int q1 = 0; q1 <= Math.min(N, maxVal); q1++) {
            for (int q2 = 0; q2 <= Math.min(N - q1, maxVal); q2++) {
                int q3 = N - q1 - q2;
                if (q3 >= 0 && q3 <= maxVal) {
                    result.add(new int[]{q1, q2, q3});
                }
            }
        }
        return result;
    }
}
