package agh.finiteelementsmethod;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class BorderConditions {
    List<Node> nodes;
    List<Element> elements;
    GlobalData globalData;
    UniversalElement universalElement;
    HMatrix hMatrix;

    Grid grid;

    public BorderConditions(List<Node> nodes, List<Element> elements, GlobalData globalData, UniversalElement universalElement, Grid grid, HMatrix hMatrix) {
        this.nodes = nodes;
        this.elements = elements;
        this.globalData = globalData;
        this.universalElement = universalElement;
        this.grid = grid;
        this.hMatrix = hMatrix;
    }

    // Gauss nodes for 2 points
    private final Point[] pointArray2P = {
            new Point(-1/sqrt(3), -1),
            new Point(1/sqrt(3), -1),
            new Point(1, -1/sqrt(3)),
            new Point(1, 1/sqrt(3)),
            new Point(-1/sqrt(3), 1),
            new Point(1/sqrt(3), 1),
            new Point(-1, -1/sqrt(3)),
            new Point(-1, 1/sqrt(3))
    };

    // Gauss nodes for 3 points
    private final Point[] pointArray3P = {
            new Point(-sqrt(3.0/5.0), -1),
            new Point(0, -1),
            new Point(sqrt(3.0/5.0), -1),
            new Point(1, -sqrt(3.0/5.0)),
            new Point(1, 0),
            new Point(1, sqrt(3.0/5.0)),
            new Point(-sqrt(3.0/5.0), 1),
            new Point(0, 1),
            new Point(sqrt(3.0/5.0), 1),
            new Point(-1, -sqrt(3.0/5.0)),
            new Point(-1, 0),
            new Point(-1, sqrt(3.0/5.0)),
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


    List<int[]> getElementsNodes() {
        List<int[]> elementsNodes = new ArrayList<>();
        for (Element element : elements) {
            int[] elementTable = new int[4];
            for (int i = 0; i < 4; i++) {
                elementTable[i] = element.getID(i);
            }
            elementsNodes.add(elementTable);
        }
        return elementsNodes;
    }

    // --------------------------------------------------------------------------------- Matrix H with border conditions

    double[][] calculateBorderConditions(int points) {
        int size = points * points;

        HMatrix hMatrix = new HMatrix(globalData, grid, nodes, elements);
        double[][] resultArray = Agregation.createGlobalMatrix(nodes);
        List<int[]> elementsNodes = getElementsNodes();
        Map<Integer, Node> nodesMap =
                nodes.stream().collect(Collectors.toMap(Node::getId, node -> node));
        int elementID = 0;

        // List of matrixes H for elements
        List<double[][]> matrixesHLocalList = new ArrayList<>();

        for (int[] element : elementsNodes) {
            double[][] elementArray1 = new double[4][4];
            double[][] elementArray2 = new double[4][4];
            double[][] elementArray3 = new double[4][4];
            double[][] elementArray4 = new double[4][4];
            double[][] elementLocalH = hMatrix.calculateMatrixH(points, elementID);

            elementID++;

            // Calculating sites for 2 points
            if (points == 2) {
                // Calculating first wall
                if (nodesMap.get(element[0]).getBC() == 1 && nodesMap.get(element[1]).getBC() == 1) {
                    elementArray1 = addTwoArrays(elementArray1,
                            calculateSide2P(0, 1, nodesMap.get(element[0]), nodesMap.get(element[1])),
                            4);
                }
                // Calculating second wall
                if (nodesMap.get(element[1]).getBC() == 1 && nodesMap.get(element[2]).getBC() == 1) {
                    elementArray2 = addTwoArrays(elementArray2,
                            calculateSide2P(2, 3, nodesMap.get(element[1]), nodesMap.get(element[2])),
                            4);
                }
                // Calculating third wall
                if (nodesMap.get(element[2]).getBC() == 1 && nodesMap.get(element[3]).getBC() == 1) {
                    elementArray3 = addTwoArrays(elementArray3,
                            calculateSide2P(4, 5, nodesMap.get(element[2]), nodesMap.get(element[3])),
                            4);
                }
                // Calculating fourth wall
                if (nodesMap.get(element[3]).getBC() == 1 && nodesMap.get(element[0]).getBC() == 1) {
                    elementArray4 = addTwoArrays(elementArray4,
                            calculateSide2P(6, 7, nodesMap.get(element[3]), nodesMap.get(element[0])),
                            4);
                }
            }

            // Calculating sites for 3 points
            if (points == 3) {
                // Calculating first wall
                if (nodesMap.get(element[0]).getBC() == 1 && nodesMap.get(element[1]).getBC() == 1) {
                    elementArray1 = addTwoArrays(elementArray1,
                            calculateSide3P(0, 1, 2, nodesMap.get(element[0]), nodesMap.get(element[1])),
                            4);
                }
                // Calculating second wall
                if (nodesMap.get(element[1]).getBC() == 1 && nodesMap.get(element[2]).getBC() == 1) {
                    elementArray2 = addTwoArrays(elementArray2,
                            calculateSide3P(3, 4, 5, nodesMap.get(element[1]), nodesMap.get(element[2])),
                            4);
                }
                // Calculating third wall
                if (nodesMap.get(element[2]).getBC() == 1 && nodesMap.get(element[3]).getBC() == 1) {
                    elementArray3 = addTwoArrays(elementArray3,
                            calculateSide3P(6, 7, 8, nodesMap.get(element[2]), nodesMap.get(element[3])),
                            4);
                }
                // Calculating fourth wall
                if (nodesMap.get(element[3]).getBC() == 1 && nodesMap.get(element[0]).getBC() == 1) {
                    elementArray4 = addTwoArrays(elementArray4,
                            calculateSide3P(9, 10, 11, nodesMap.get(element[3]), nodesMap.get(element[0])),
                            4);
                }
            }

            double[][] addedSidesMatrix = addFourArrays(elementArray1, elementArray2, elementArray3, elementArray4, 4);
            double[][] elementHResult = addTwoArrays(elementLocalH, addedSidesMatrix, 4);
            matrixesHLocalList.add(elementHResult);
        }
        assignLocalHToElements(matrixesHLocalList);
        double[][] temporaryGlobalH = Agregation.createGlobalMatrix(nodes);
        resultArray = Agregation.agregateH(temporaryGlobalH, elements);
        return resultArray;
    }

    private void assignLocalHToElements(List<double[][]> localHList) {
        int counter = 0;
        for (Element element : elements) {
            element.setLocalH(localHList.get(counter));
            counter++;
        }
    }

    double[][] addTwoArrays(double[][] firstArray, double[][] secondArray, int size) {
        double[][] resultArray = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                resultArray[i][j] += firstArray[i][j];
                resultArray[i][j] += secondArray[i][j];
            }
        }
        return resultArray;
    }

    double[][] addFourArrays(double[][] array1, double[][] array2, double[][] array3, double[][] array4, int size) {
        double[][] resultArray = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                resultArray[i][j] += array1[i][j];
                resultArray[i][j] += array2[i][j];
                resultArray[i][j] += array3[i][j];
                resultArray[i][j] += array4[i][j];
            }
        }
        return resultArray;
    }

    double[][] calculateSide2P(int firstPoint, int secondPoint, Node firstNode, Node secondNode) {
        double[][] side = new double[4][4];
        double[] shapeFunctionValuesFirstPoint = new double[] {
          universalElement.N1(pointArray2P[firstPoint].getKsi(), pointArray2P[firstPoint].getEta()),
          universalElement.N2(pointArray2P[firstPoint].getKsi(), pointArray2P[firstPoint].getEta()),
          universalElement.N3(pointArray2P[firstPoint].getKsi(), pointArray2P[firstPoint].getEta()),
          universalElement.N4(pointArray2P[firstPoint].getKsi(), pointArray2P[firstPoint].getEta())
        };
        double[] shapeFunctionValuesSecondPoint = new double[] {
            universalElement.N1(pointArray2P[secondPoint].getKsi(), pointArray2P[secondPoint].getEta()),
            universalElement.N2(pointArray2P[secondPoint].getKsi(), pointArray2P[secondPoint].getEta()),
            universalElement.N3(pointArray2P[secondPoint].getKsi(), pointArray2P[secondPoint].getEta()),
            universalElement.N4(pointArray2P[secondPoint].getKsi(), pointArray2P[secondPoint].getEta())
        };
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                // Wages are skipped because they are equal to 1
                side[i][j] = globalData.getAlfa() * (shapeFunctionValuesFirstPoint[i] * shapeFunctionValuesFirstPoint[j] + shapeFunctionValuesSecondPoint[i] * shapeFunctionValuesSecondPoint[j])
                        * calculateDetJ(firstNode, secondNode);
            }
        }

        return side;
    }

    double[][] calculateSide3P(int firstPoint, int secondPoint, int thirdPoint, Node firstNode, Node secondNode) {
        double[][] side = new double[4][4];
        double[] shapeFunctionValuesFirstPoint = new double[] {
                universalElement.N1(pointArray3P[firstPoint].getKsi(), pointArray3P[firstPoint].getEta()),
                universalElement.N2(pointArray3P[firstPoint].getKsi(), pointArray3P[firstPoint].getEta()),
                universalElement.N3(pointArray3P[firstPoint].getKsi(), pointArray3P[firstPoint].getEta()),
                universalElement.N4(pointArray3P[firstPoint].getKsi(), pointArray3P[firstPoint].getEta())
        };
        double[] shapeFunctionValuesSecondPoint = new double[] {
                universalElement.N1(pointArray3P[secondPoint].getKsi(), pointArray3P[secondPoint].getEta()),
                universalElement.N2(pointArray3P[secondPoint].getKsi(), pointArray3P[secondPoint].getEta()),
                universalElement.N3(pointArray3P[secondPoint].getKsi(), pointArray3P[secondPoint].getEta()),
                universalElement.N4(pointArray3P[secondPoint].getKsi(), pointArray3P[secondPoint].getEta())
        };
        double[] shapeFunctionValuesThirdPoint = new double[] {
                universalElement.N1(pointArray3P[thirdPoint].getKsi(), pointArray3P[thirdPoint].getEta()),
                universalElement.N2(pointArray3P[thirdPoint].getKsi(), pointArray3P[thirdPoint].getEta()),
                universalElement.N3(pointArray3P[thirdPoint].getKsi(), pointArray3P[thirdPoint].getEta()),
                universalElement.N4(pointArray3P[thirdPoint].getKsi(), pointArray3P[thirdPoint].getEta())
        };
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                side[i][j] = globalData.getAlfa() *
                        (shapeFunctionValuesFirstPoint[i] * shapeFunctionValuesFirstPoint[j] * 5/9
                                + shapeFunctionValuesSecondPoint[i] * shapeFunctionValuesSecondPoint[j] * 8/9
                        + shapeFunctionValuesThirdPoint[i] * shapeFunctionValuesThirdPoint[j] * 5/9
                        )
                        * calculateDetJ(firstNode, secondNode);
            }
        }

        return side;
    }

    double calculateDetJ(Node node1, Node node2) {
        return (sqrt(pow(node1.getX() - node2.getX(), 2) + pow(node1.getY() - node2.getY(), 2))) / 2;
    }

    // -----------------------------------------------------------------------------------------------------------------

    // -------------------------------------------------------------------------------------------------------- Vector P

    double[] calculateVectorPPart2P(int firstPoint, int secondPoint, Node firstNode, Node secondNode) {
        double[] vectorPPart = new double[4];
        double[] shapeFunctionValuesFirstPoint = new double[] {
                universalElement.N1(pointArray2P[firstPoint].getKsi(), pointArray2P[firstPoint].getEta()),
                universalElement.N2(pointArray2P[firstPoint].getKsi(), pointArray2P[firstPoint].getEta()),
                universalElement.N3(pointArray2P[firstPoint].getKsi(), pointArray2P[firstPoint].getEta()),
                universalElement.N4(pointArray2P[firstPoint].getKsi(), pointArray2P[firstPoint].getEta())
        };
        double[] shapeFunctionValuesSecondPoint = new double[] {
                universalElement.N1(pointArray2P[secondPoint].getKsi(), pointArray2P[secondPoint].getEta()),
                universalElement.N2(pointArray2P[secondPoint].getKsi(), pointArray2P[secondPoint].getEta()),
                universalElement.N3(pointArray2P[secondPoint].getKsi(), pointArray2P[secondPoint].getEta()),
                universalElement.N4(pointArray2P[secondPoint].getKsi(), pointArray2P[secondPoint].getEta())
        };

        for (int i = 0; i < 4; i++) {
            // Wages are skipped because they are equal to 1
            vectorPPart[i] = globalData.getAlfa() * (shapeFunctionValuesFirstPoint[i] * globalData.getTot() + shapeFunctionValuesSecondPoint[i] * globalData.getTot())
                    * calculateDetJ(firstNode, secondNode);
        }
        return vectorPPart;
    }

    double[] calculateVectorPPart3P(int firstPoint, int secondPoint, int thirdPoint, Node firstNode, Node secondNode) {
        double[] vectorPPart = new double[4];
        double[] shapeFunctionValuesFirstPoint = new double[] {
                universalElement.N1(pointArray3P[firstPoint].getKsi(), pointArray3P[firstPoint].getEta()),
                universalElement.N2(pointArray3P[firstPoint].getKsi(), pointArray3P[firstPoint].getEta()),
                universalElement.N3(pointArray3P[firstPoint].getKsi(), pointArray3P[firstPoint].getEta()),
                universalElement.N4(pointArray3P[firstPoint].getKsi(), pointArray3P[firstPoint].getEta())
        };
        double[] shapeFunctionValuesSecondPoint = new double[] {
                universalElement.N1(pointArray3P[secondPoint].getKsi(), pointArray3P[secondPoint].getEta()),
                universalElement.N2(pointArray3P[secondPoint].getKsi(), pointArray3P[secondPoint].getEta()),
                universalElement.N3(pointArray3P[secondPoint].getKsi(), pointArray3P[secondPoint].getEta()),
                universalElement.N4(pointArray3P[secondPoint].getKsi(), pointArray3P[secondPoint].getEta())
        };
        double[] shapeFunctionValuesThirdPoint = new double[] {
                universalElement.N1(pointArray3P[thirdPoint].getKsi(), pointArray3P[thirdPoint].getEta()),
                universalElement.N2(pointArray3P[thirdPoint].getKsi(), pointArray3P[thirdPoint].getEta()),
                universalElement.N3(pointArray3P[thirdPoint].getKsi(), pointArray3P[thirdPoint].getEta()),
                universalElement.N4(pointArray3P[thirdPoint].getKsi(), pointArray3P[thirdPoint].getEta())
        };

        for (int i = 0; i < 4; i++) {
            // Wages are skipped because they are equal to 1
            vectorPPart[i] = globalData.getAlfa() *
                    (shapeFunctionValuesFirstPoint[i] * globalData.getTot() * 5/9
                            + shapeFunctionValuesSecondPoint[i] * globalData.getTot() * 8/9
                    + shapeFunctionValuesThirdPoint[i] * globalData.getTot() * 5/9
                    )
                    * calculateDetJ(firstNode, secondNode);
        }
        return vectorPPart;
    }

    double[] calculateVectorP(int points) {
        double[] globalVectorP = Agregation.createGlobalVectorP(nodes);
        List<int[]> elementsNodes = getElementsNodes();
        Map<Integer, Node> nodesMap =
                nodes.stream().collect(Collectors.toMap(Node::getId, node -> node));

        List<double[]> vectorsP = new ArrayList<>();

        for (int[] element : elementsNodes) {
            double[] part1VectorP = new double[4];
            double[] part2VectorP = new double[4];
            double[] part3VectorP = new double[4];
            double[] part4VectorP = new double[4];

            if (points == 2) {
                // Calculating first wall
                if (nodesMap.get(element[0]).getBC() == 1 && nodesMap.get(element[1]).getBC() == 1) {
                    part1VectorP = addTwoVectors(part1VectorP,
                            calculateVectorPPart2P(0, 1, nodesMap.get(element[0]), nodesMap.get(element[1])),
                            4);
                }
                // Calculating second wall
                if (nodesMap.get(element[1]).getBC() == 1 && nodesMap.get(element[2]).getBC() == 1) {
                    part2VectorP = addTwoVectors(part2VectorP,
                            calculateVectorPPart2P(2, 3, nodesMap.get(element[1]), nodesMap.get(element[2])),
                            4);
                }
                // Calculating third wall
                if (nodesMap.get(element[2]).getBC() == 1 && nodesMap.get(element[3]).getBC() == 1) {
                    part3VectorP = addTwoVectors(part3VectorP,
                            calculateVectorPPart2P(4, 5, nodesMap.get(element[2]), nodesMap.get(element[3])),
                            4);
                }
                // Calculating fourth wall
                if (nodesMap.get(element[3]).getBC() == 1 && nodesMap.get(element[0]).getBC() == 1) {
                    part4VectorP = addTwoVectors(part4VectorP,
                            calculateVectorPPart2P(6, 7, nodesMap.get(element[3]), nodesMap.get(element[0])),
                            4);
                }
            }

            if (points == 3) {
                // Calculating first wall
                if (nodesMap.get(element[0]).getBC() == 1 && nodesMap.get(element[1]).getBC() == 1) {
                    part1VectorP = addTwoVectors(part1VectorP,
                            calculateVectorPPart3P(0, 1, 2, nodesMap.get(element[0]), nodesMap.get(element[1])),
                            4);
                }
                // Calculating second wall
                if (nodesMap.get(element[1]).getBC() == 1 && nodesMap.get(element[2]).getBC() == 1) {
                    part2VectorP = addTwoVectors(part2VectorP,
                            calculateVectorPPart3P(3, 4, 5, nodesMap.get(element[1]), nodesMap.get(element[2])),
                            4);
                }
                // Calculating third wall
                if (nodesMap.get(element[2]).getBC() == 1 && nodesMap.get(element[3]).getBC() == 1) {
                    part3VectorP = addTwoVectors(part3VectorP,
                            calculateVectorPPart3P(6, 7, 8, nodesMap.get(element[2]), nodesMap.get(element[3])),
                            4);
                }
                // Calculating fourth wall
                if (nodesMap.get(element[3]).getBC() == 1 && nodesMap.get(element[0]).getBC() == 1) {
                    part4VectorP = addTwoVectors(part4VectorP,
                            calculateVectorPPart3P(9, 10, 11, nodesMap.get(element[3]), nodesMap.get(element[0])),
                            4);
                }
            }

            double[] localVectorP = addFourVectors(part1VectorP, part2VectorP, part3VectorP, part4VectorP, 4);
            vectorsP.add(localVectorP);
        }

        double[] temporaryGlobalVectorP = Agregation.createGlobalVectorP(nodes);
        globalVectorP = Agregation.agregateVectorP(vectorsP, temporaryGlobalVectorP, elements);

        return globalVectorP;
    }

    double[] addTwoVectors(double[] firstVector, double[] secondVector, int size) {
        double[] resultVector = new double[size];
        for (int i = 0; i < size; i++) {
            resultVector[i] += firstVector[i];
            resultVector[i] += secondVector[i];
        }
        return resultVector;
    }

    double[] addFourVectors(double[] vector1, double[] vector2, double[] vector3, double[] vector4, int size) {
        double[] resultVector = new double[size];
        for (int i = 0; i < size; i++) {
            resultVector[i] += vector1[i];
            resultVector[i] += vector2[i];
            resultVector[i] += vector3[i];
            resultVector[i] += vector4[i];
        }
        return resultVector;
    }

    // -----------------------------------------------------------------------------------------------------------------

    // -------------------------------------------------------------------------------------------------------- Matrix C

    double[][] calculateMatrixC(int points) {
        CMatrix cMatrix = new CMatrix(globalData, grid, nodes, elements);
        double[][] resultArray = Agregation.createGlobalMatrix(nodes);
        List<int[]> elementsNodes = getElementsNodes();
        Map<Integer, Node> nodesMap =
                nodes.stream().collect(Collectors.toMap(Node::getId, node -> node));
        int elementID = 0;

        // List of matrixes C for elements
        List<double[][]> matrixesCLocalList = new ArrayList<>();

        for (int[] element : elementsNodes) {
            double[][] elementLocalC = cMatrix.calculateMatrixC(points, elementID);
            elementID++;

            // Adding every local C to a list
            matrixesCLocalList.add(elementLocalC);
        }
        assignLocalCToElements(matrixesCLocalList);
        double[][] temporaryGlobalC = Agregation.createGlobalMatrix(nodes);
        resultArray = Agregation.agregateC(temporaryGlobalC, elements);
        return resultArray;
    }

    private void assignLocalCToElements(List<double[][]> localCList) {

        int counter = 0;
        for (Element element : elements) {
            element.setLocalC(localCList.get(counter));
            counter++;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------

}
