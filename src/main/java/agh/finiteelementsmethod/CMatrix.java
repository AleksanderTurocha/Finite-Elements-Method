package agh.finiteelementsmethod;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;

import static java.lang.Math.sqrt;

public class CMatrix {

    // Global data is used here to pass conductivity during calculating final stage of H matrix
    private final GlobalData globalData;
    private final Grid grid;

    private final List<Node> nodes;

    private final List<Element> elements;

    public CMatrix(GlobalData globalData, Grid grid, List<Node> nodes, List<Element> elements)
    {
        this.globalData = globalData;
        this.grid = grid;
        this.nodes = nodes;
        this.elements = elements;
    }

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
    double[] wages2P = {1, 1, 1, 1};

    //Wages for 3 points
    double[] wages3P = {5.0/9.0, 8.0/9.0, 5.0/9.0};

    //Wages for 4 points
    double[] wages4P = {(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36};

    public static double[] exportDetJArray(double[] detJArray) {
        return detJArray;
    }

    // Main calculating function
    double[][] calculateMatrixC(int points, int elementID) {

        // Element nodes dropped into the list
        Element wantedElement = elements.get(elementID);
        Node elementNode1 = nodes.get(wantedElement.getID(0)-1);
        Node elementNode2 = nodes.get(wantedElement.getID(1)-1);
        Node elementNode3 = nodes.get(wantedElement.getID(2)-1);
        Node elementNode4 = nodes.get(wantedElement.getID(3)-1);

        List<Node> elementNodesList = List.of(elementNode1, elementNode2, elementNode3, elementNode4);

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
        double[][] arrayJ = calculateJ(size, universalElement, elementNodesList);

        // Calculating detJ
        double[] detJArray = new double[size];
        double[] detJArrayReversed = new double[size];

        for (int i = 0; i < size; i++) {
            detJArray[i] = arrayJ[i][0] * arrayJ[i][3] - arrayJ[i][1] * arrayJ[i][2];
        }

        for (int i = 0; i < size; i++) {
            detJArrayReversed[i] = 1/detJArray[i];
        }

        // Calculating arrayC - array of matrixes C
        double[][][] arrayC = calculateArrayC(size, detJArray, pointArray);

        // Applying wages
        Double[] wagesArray = getWages(wages, points);

        // Applying wages for matrix C
        double[][] C = applyWages(arrayC, size, wagesArray);
        return C;
    }

    double[][][] calculateArrayC(int size, double[] detJArray, Point[] pointArray) {
        double[][][] arrayC = new double[size][4][4];
        List<BiFunction<Double, Double, Double>> shapeFunctions = new ArrayList<>() {{
            add((ksi, eta) -> 0.25 * (1 - ksi) * (1 - eta));
            add((ksi, eta) -> 0.25 * (1 + ksi) * (1 - eta));
            add((ksi, eta) -> 0.25 * (1 + ksi) * (1 + eta));
            add((ksi, eta) -> 0.25 * (1 - ksi) * (1 + eta));
        }};

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    arrayC[i][j][k] = globalData.getDensity() * globalData.getSpecificHeat() *
                shapeFunctions.get(k).apply(pointArray[i].getKsi(), pointArray[i].getEta())
                            * shapeFunctions.get(j).apply(pointArray[i].getKsi(), pointArray[i].getEta())
                * detJArray[i];
                }
            }
        }
        return arrayC;
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

    private double[][] applyWages(double[][][] arrayC , int size, Double[] wages) {
        double[][] C = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < size; k++) {
                    C[i][j] += arrayC[k][i][j] * wages[k];
                }
            }
        }
        return C;
    }

    void displayC(double[][] C) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                System.out.printf("|%-6.2f|", C[j][i]);
            }
            System.out.println();
        }
    }

    private double[][] calculateJ(int size, UniversalElement universalElement, List<Node> elementNodesList) {
        double[][] arrayJ = new double[size][4];
        for (int j = 0; j < size; j++) {
            double dxdKsi = 0;
            double dydKsi = 0;
            double dxdEta = 0;
            double dydEta = 0;
            for (int i = 0; i < 4; i++) {
                dxdKsi += elementNodesList.get(i).getX() * universalElement.getArraydNdKsi()[i][j];
                dydKsi += elementNodesList.get(i).getY() * universalElement.getArraydNdKsi()[i][j];
                dxdEta += elementNodesList.get(i).getX() * universalElement.getArraydNdEta()[i][j];
                dydEta += elementNodesList.get(i).getY() * universalElement.getArraydNdEta()[i][j];
                arrayJ[j][0] = dxdKsi;
                arrayJ[j][1] = dydKsi;
                arrayJ[j][2] = dxdEta;
                arrayJ[j][3] = dydEta;
            }
        }
        return arrayJ;
    }

    private double[][] calculateDNDEta(int size, Point[] pointArray) {
        UniversalElement universalElement = new UniversalElement();
        double[][] arraydNdEta = new double[4][size];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < size; j++) {
                arraydNdEta[i][j] = universalElement.dNdEtaFunctions(i, pointArray[j].getKsi());
            }
        }
        return arraydNdEta;
    }

    private double[][] calculateDNDKsi(int size, Point[] pointArray) {
        UniversalElement universalElement = new UniversalElement();
        double[][] arraydNdKsi = new double[4][size];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < size; j++) {
                arraydNdKsi[i][j] = universalElement.dNdKsiFunctions(i, pointArray[j].getEta());
            }
        }
        return arraydNdKsi;
    }

}

