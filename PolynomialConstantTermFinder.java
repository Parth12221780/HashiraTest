import java.io.*;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import org.json.*;

public class PolynomialConstantTermFinder {

    // Convert encoded string with given base to BigInteger
    public static BigInteger decodeValue(String value, int base) {
        return new BigInteger(value, base);
    }

    // Lagrange interpolation to find f(0)
    // points is list of (x, y) pairs where x,y are BigInteger
    public static BigInteger lagrangeInterpolationAtZero(List<BigInteger[]> points) {
        BigInteger result = BigInteger.ZERO;

        for (int j = 0; j < points.size(); j++) {
            BigInteger xj = points.get(j)[0];
            BigInteger yj = points.get(j)[1];

            BigInteger numerator = BigInteger.ONE;
            BigInteger denominator = BigInteger.ONE;

            for (int m = 0; m < points.size(); m++) {
                if (m != j) {
                    BigInteger xm = points.get(m)[0];
                    numerator = numerator.multiply(xm).mod(null); // no mod, keep normal multiplication
                    denominator = denominator.multiply(xm.subtract(xj)).mod(null);
                }
            }
            
            BigInteger numeratorProduct = yj.multiply(numerator);
            // Divide by denominator (BigInteger division)
            BigInteger term = numeratorProduct.divide(denominator);
            result = result.add(term);
        }

        return result;
    }

   

    // Solve system of linear equations A * coeffs = Y using Gaussian elimination
    // A is k x k matrix, coeffs is k x 1 vector of coefficients [a0=c, a1,...,a_m]
    // Y is k x 1 vector of y values
    public static BigInteger[] solveCoefficients(List<BigInteger[]> points, int k) {
        // Build matrix A and vector Y
        BigInteger[][] A = new BigInteger[k][k];
        BigInteger[] Y = new BigInteger[k];

        for (int i = 0; i < k; i++) {
            BigInteger x = points.get(i)[0];
            Y[i] = points.get(i)[1];
            BigInteger pow = BigInteger.ONE;
            for (int j = 0; j < k; j++) {
                A[i][j] = pow;
                pow = pow.multiply(x);
            }
        }

        // Solve linear system A * coeffs = Y
        return gaussianElimination(A, Y);
    }

    // Gaussian elimination for BigInteger matrix
    public static BigInteger[] gaussianElimination(BigInteger[][] A, BigInteger[] Y) {
        int n = Y.length;

        for (int i = 0; i < n; i++) {
            // Search for maximum in this column
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (A[k][i].abs().compareTo(A[maxRow][i].abs()) > 0) {
                    maxRow = k;
                }
            }

            // Swap max row with current row
            BigInteger[] tempRow = A[i];
            A[i] = A[maxRow];
            A[maxRow] = tempRow;

            BigInteger tempVal = Y[i];
            Y[i] = Y[maxRow];
            Y[maxRow] = tempVal;

            // Make all rows below this one 0 in current column
            for (int k = i + 1; k < n; k++) {
                if (A[k][i].equals(BigInteger.ZERO)) continue;

                // Calculate factor = A[k][i] / A[i][i]
                // Division of BigInteger not exact; we'll do fraction carefully
                BigInteger numerator = A[k][i];
                BigInteger denominator = A[i][i];

                // We can't do fractional arithmetic directly,
                // so we multiply entire row i by numerator and row k by denominator and subtract
                for (int j = i; j < n; j++) {
                    A[k][j] = A[k][j].multiply(denominator).subtract(A[i][j].multiply(numerator));
                }
                Y[k] = Y[k].multiply(denominator).subtract(Y[i].multiply(numerator));
            }
        }

        // Back substitution
        BigInteger[] x = new BigInteger[n];
        for (int i = n - 1; i >= 0; i--) {
            BigInteger sum = BigInteger.ZERO;
            for (int j = i + 1; j < n; j++) {
                sum = sum.add(A[i][j].multiply(x[j]));
            }
            // x[i] = (Y[i] - sum) / A[i][i]
            x[i] = Y[i].subtract(sum).divide(A[i][i]);
        }
        return x;
    }

    public static void processFile(String filePath) throws IOException {
        String content = new String(Files.readAllBytes(Paths.get(filePath)));
        JSONObject json = new JSONObject(content);

        JSONObject keys = json.getJSONObject("keys");
        int n = keys.getInt("n");
        int k = keys.getInt("k");

        // Parse points
        List<BigInteger[]> points = new ArrayList<>();

        for (String key : json.keySet()) {
            if (key.equals("keys")) continue;
            int x = Integer.parseInt(key);
            JSONObject rootObj = json.getJSONObject(key);

            int base = Integer.parseInt(rootObj.getString("base"));
            String value = rootObj.getString("value");

            BigInteger y = new BigInteger(value, base);
            points.add(new BigInteger[]{BigInteger.valueOf(x), y});
        }

        // Sort by x
        points.sort(Comparator.comparing(p -> p[0]));

        // Use first k points only
        List<BigInteger[]> usedPoints = points.subList(0, k);

        BigInteger[] coeffs = solveCoefficients(usedPoints, k);

        // constant term c = coeffs[0]
        System.out.println(coeffs[0]);
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.out.println("Usage: java PolynomialConstantTermFinder <file1.json> <file2.json>");
            return;
        }
        processFile(args[0]);
        processFile(args[1]);
    }
}
