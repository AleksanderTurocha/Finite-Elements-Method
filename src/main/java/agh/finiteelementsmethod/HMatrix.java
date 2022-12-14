package agh.finiteelementsmethod;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.round;
import static java.lang.Math.sqrt;

public class HMatrix {
    // Temperature set at the beginning
    private static final double CONDUCT = 30;

    // Gauss nodes for 2 points
    private final Point[] pointArray2P = {
            new Point(-1/sqrt(3), -1/sqrt(3)),
            new Point(1/sqrt(3), -1/sqrt(3)),
            new Point(-1/sqrt(3), 1/sqrt(3)),
            new Point(1/sqrt(3), 1/sqrt(3))
    };

    // Gauss nodes for 3 points
    private final Point[] pointArray3P = {
            new Point(-sqrt(3.0/5.0), -sqrt(3.0/5.0)),
            new Point(0, -sqrt(3.0/5.0)),
            new Point(sqrt(3.0/5.0), -sqrt(3.0/5.0)),
            new Point(-sqrt(3.0/5.0), 0),
            new Point(0, 0),
            new Point(sqrt(3.0/5.0), 0),
            new Point(-sqrt(3.0/5.0), sqrt(3.0/5.0)),
            new Point(0, sqrt(3.0/5.0)),
            new Point(sqrt(3.0/5.0), sqrt(3.0/5.0))
    };

    //  Gauss nodes for 4 points
    private final Point[] pointArray4P = {
        new Point (-sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5)), -sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5))),
        new Point(-sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5)), -sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5))),
        new Point(sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5)), -sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5))),
        new Point(sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5)), -sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5))),
        new Point(-sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5)), -sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5))),
        new Point(-sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5)), -sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5))),
        new Point(sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5)), -sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5))),
        new Point(sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5)), -sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5))),
        new Point(-sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5)), sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5))),
        new Point(-sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5)), sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5))),
        new Point(sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5)), sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5))),
        new Point(sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5)), sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5))),
        new Point(-sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5)), sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5))),
        new Point(-sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5)), sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5))),
        new Point(sqrt( 3.0/7 - 2.0/7 * sqrt(6.0/5)), sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5))),
        new Point(sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5)), sqrt( 3.0/7 + 2.0/7 * sqrt(6.0/5.0)))
    };

    // Wages for 2 points
    double[] wages2P = {1, 1};

    //Wages for 3 points
    double[] wages3P = {5.0/9.0, 8.0/9.0, 5.0/9.0};

    //Wages for 4 points
    double[] wages4P = {(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36};

    // Point subclass holding Gauss node coordinates
    static class Point {
        private double ksi;
        private double eta;

        public double getKsi() {
            return ksi;
        }

        public void setKsi(double ksi) {
            this.ksi = ksi;
        }

        public double getEta() {
            return eta;
        }

        public void setEta(double eta) {
            this.eta = eta;
        }

        public Point(double ksi, double eta) {
            this.ksi = ksi;
            this.eta = eta;
        }
    }

    // Counting derivatives of shape functions
    // Interface of a derivative a shape function
    public interface DShapeFunction {
        double calculate(double value);
    }

    // Implementation of all derivative functions in form of a class with an interface
    // N1 / ksi
    static class DN1DKsi implements DShapeFunction {
        @Override
        public double calculate(double eta) {
            return -0.25 * (1 - eta);
        }
    }

    // N2 / ksi
    static class DN2DKsi implements DShapeFunction {
        @Override
        public double calculate(double eta) {
            return 0.25 * (1 - eta);
        }
    }

    // N3 / ksi
    static class DN3DKsi implements DShapeFunction {
        @Override
        public double calculate(double eta) {
            return 0.25 * (1 + eta);
        }
    }

    // N4 / ksi
    static class DN4DKsi implements DShapeFunction {
        @Override
        public double calculate(double eta) {
            return -0.25 * (1 + eta);
        }
    }

    // N1 / eta
    static class DN1DEta implements DShapeFunction {
        @Override
        public double calculate(double ksi) {
            return -0.25 * (1 - ksi);
        }
    }

    // N2 / eta
    static class DN2DEta implements DShapeFunction {
        @Override
        public double calculate(double ksi) {
            return -0.25 * (1 + ksi);
        }
    }

    // N3 / eta
    static class DN3DEta implements DShapeFunction {
        @Override
        public double calculate(double ksi) {
            return 0.25 * (1 + ksi);
        }
    }

    // N4 / eta
    static class DN4DEta implements DShapeFunction {
        @Override
        public double calculate(double ksi) {
            return 0.25 * (1 - ksi);
        }
    }

    // Creating tables derivatives eta and ksi
    static final DShapeFunction[] dNdKsiFunctions = {new DN1DKsi(), new DN2DKsi(), new DN3DKsi(), new DN4DKsi()};
    static final DShapeFunction[] dNdEtaFunctions = {new DN1DEta(), new DN2DEta(), new DN3DEta(), new DN4DEta()};

    // Universal Element structure containing arrays of dNdKsi and dNdEta
    static class UniversalElement {
        private double[][] arraydNdKsi;
        private double[][] arraydNdEta;

        public UniversalElement(double[][] arraydNdKsi, double[][] arraydNdEta) {
            this.arraydNdKsi = arraydNdKsi;
            this.arraydNdEta = arraydNdEta;
        }


        public double[][] getArraydNdKsi() {
            return arraydNdKsi;
        }

        public void setArraydNdKsi(double[][] arraydNdKsi) {
            this.arraydNdKsi = arraydNdKsi;
        }

        public double[][] getArraydNdEta() {
            return arraydNdEta;
        }

        public void setArraydNdEta(double[][] arraydNdEta) {
            this.arraydNdEta = arraydNdEta;
        }
    }

    // Given point coordinates dropped into the list
    static Point point1 = new Point(0, 0);
    static Point point2 = new Point(0.025, 0);
    static Point point3 = new Point(0.025, 0.025);
    static Point point4 = new Point(0, 0.025);

    static List<Point> pointList = List.of(point1, point2, point3, point4);

    // Main calculating function
    double[][] calculatingMatrixH (int points) {
        int size = points * points;
        double[][] arraydNdKsi = new double[4][size];
        double[][] arraydNdEta = new double[4][size];
        Point[] pointArray = new Point[size];
        double[] wages = new double[size];
        if (points == 2) {
            pointArray = pointArray2P;
            wages = wages2P;
        } else if (points == 3) {
            pointArray = pointArray3P;
            wages = wages3P;
        } else if (points == 4) {
            pointArray = pointArray4P;
            wages = wages4P;
        }

        // Calculating dNdKsi
        arraydNdKsi = calculateDNDKsi(size, pointArray);

        // Calculating dNdEta
        arraydNdEta = calculateDNDEta(size, pointArray);

        // Parsing actual dNdEta and dNdKsi arrays to the universal element object
        UniversalElement universalElement = new UniversalElement(arraydNdKsi, arraydNdEta);

        // Array J = dydEta and dxdKsi
        double[][] arrayJ = calculateJ(size, universalElement);

        // Calculating detJ
        double[] detJArray = new double[size];
        double[] detJArrayReversed = new double[size];

        for (int i = 0; i < size; i++) {
            detJArray[i] = arrayJ[i][0] * arrayJ[i][3] - arrayJ[i][1] * arrayJ[i][2];
        }

        for (int i = 0; i < size; i++) {
            detJArrayReversed[i] = 1/detJArray[i];
        }

        // Derivatives of x and y functions / ksi and eta
        double[][][] dNdxArray = new double[size][4][size];
        double[][][] dNdyArray = new double[size][4][size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < size; k++) {
                    dNdxArray[i][j][k] = detJArrayReversed[i] * arrayJ[i][3] * universalElement.getArraydNdKsi()[j][k] + detJArrayReversed[i] * (-arrayJ[i][1]) * universalElement.getArraydNdEta()[j][k];
                }
                for (int l = 0; l < size; l++) {
                    dNdyArray[i][j][l] = detJArrayReversed[i] * (-arrayJ[i][2]) * universalElement.getArraydNdKsi()[j][l] + detJArrayReversed[i] * arrayJ[i][0] * universalElement.getArraydNdEta()[j][l];
                }
            }
        }

        // Creating H1, H2 etc. (parts of H matrix) - multiplying values from tables dNdx and dNdy like vectors
        double[][][] Hdx = new double[size][4][4];
        double[][][] Hdy = new double[size][4][4];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    Hdx[i][j][k] = dNdxArray[i][j][i] * dNdxArray[i][k][i];
                }
            }
        }

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    Hdy[i][j][k] = dNdyArray[i][j][i] * dNdyArray[i][k][i];
                }
            }
        }

        // Creating array of those H1, H2 etc. matrixes and multiplying values with temperature and detJ
        double[][][] arrayH = new double[size][4][4];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    arrayH[i][j][k] = CONDUCT * (Hdx[i][j][k] + Hdy[i][j][k]) * detJArray[i];
                }
            }
        }

        // Applying wages
        Double[] wagesArray = getWages(wages, points);

        // Calculating final H - our result and applying wages
        double[][] H = createH(arrayH, size, wagesArray);
        return H;
    }

    private Double[] getWages(double[] wages, int points) {
        List<Double> wagesList = new ArrayList<>();
        for (int i = 0; i < points; i++) {
            for (int j = 0; j < points; j++) {
                wagesList.add(wages[i] * wages[j]);
            }
        }
        Double[] wagesArray = new Double[wagesList.size()];
        wagesArray = wagesList.toArray(wagesArray);
        return wagesArray;
    }

    private double[][] createH(double[][][] arrayH , int size, Double[] wages) {
        double[][] H = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < size; k++) {
                    H[i][j] += arrayH[k][i][j] * wages[k];
                }
            }
        }
        return H;
    }

    void displayH(double[][] H) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.print(round(H[j][i]) + " | ");
            }
            System.out.println();
        }
    }

    private double[][] calculateJ(int size, UniversalElement universalElement) {
        double[][] arrayJ = new double[size][4];
        for (int j = 0; j < size; j++) {
            double dxdKsi = 0;
            double dydKsi = 0;
            double dxdEta = 0;
            double dydEta = 0;
            for (int i = 0; i < 4; i++) {
                dxdKsi += pointList.get(i).getKsi() * universalElement.getArraydNdKsi()[i][j];
                dydKsi += pointList.get(i).getEta() * universalElement.getArraydNdKsi()[i][j];
                dxdEta += pointList.get(i).getKsi() * universalElement.getArraydNdEta()[i][j];
                dydEta += pointList.get(i).getEta() * universalElement.getArraydNdEta()[i][j];
                arrayJ[j][0] = dxdKsi;
                arrayJ[j][1] = dydKsi;
                arrayJ[j][2] = dxdEta;
                arrayJ[j][3] = dydEta;
            }
        }
        return arrayJ;
    }

    private double[][] calculateDNDEta(int size, Point[] pointArray) {
        double[][] arraydNdEta = new double[4][size];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < size; j++) {
                arraydNdEta[i][j] = dNdEtaFunctions[i].calculate(pointArray[j].ksi);
            }
        }
        return arraydNdEta;
    }

    private double[][] calculateDNDKsi(int size, Point[] pointArray) {
        double[][] arraydNdKsi = new double[4][size];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < size; j++) {
                arraydNdKsi[i][j] = dNdKsiFunctions[i].calculate(pointArray[j].eta);
            }
        }
        return arraydNdKsi;
    }

}
