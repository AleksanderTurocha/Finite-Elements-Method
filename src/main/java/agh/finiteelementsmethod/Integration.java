package agh.finiteelementsmethod;

public class Integration {
    private final double[] pc1 = {0.57, -0.57};
    private final double[] pc2 = {0.77, -0.77, 0};
    private final double[] wages1 = {1};
    private final double[] wages2 = {0.55, 0.88, 0.55};

    double firstEquation(double value) {
        return (2 * (value * value)) + (3 * value) + (-8);
    }

    double secondEquation(double firstValue, double secondValue) {
        return (-2 * (firstValue * firstValue) * secondValue) + ((2 * firstValue) * (secondValue * secondValue)) + 10;
    }

    double integrate1D(int points) {
        double result = 0;

        if (points == 2) {
            for (int i = 0; i < 2; i++) {
                result += firstEquation(pc1[i]) * wages1[0];
            }
        } else if (points == 3) {
            for (int i = 0; i < 3; i++) {
                result += firstEquation(pc2[i]) * wages2[i];
            }
        } else {
            System.out.println("Wrong points number!");
            return 0;
        }
        return result;
    }

    double integrate2D(int points) {
        double result = 0;
        if (points == 2) {
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    result += secondEquation(pc1[i], pc1[j]) * wages1[0] * wages1[0];
                }
            }
        } else if (points == 3) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    result += secondEquation(pc2[i], pc2[j]) * wages2[i] * wages2[j];
                }
            }
        } else {
            System.out.println("Wrong points number!");
            return 0;
        }
        return result;
    }

    public static void main(String[] args) {
        Integration integration = new Integration();
        System.out.println("Integration in 1D, 2 points: " + integration.integrate1D(2));
        System.out.println("Integration in 1D, 3 points: " + integration.integrate1D(3));
        System.out.println("Integration in 2D, 2 points: " + integration.integrate2D(2));
        System.out.println("Integration in 2D, 3 points: " + integration.integrate2D(3));
    }
}
