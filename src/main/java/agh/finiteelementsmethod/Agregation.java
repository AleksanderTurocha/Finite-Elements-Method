package agh.finiteelementsmethod;

import java.util.List;

public class Agregation {
    static double[][] createGlobalMatrix(List<Node> nodes) {
        int nodesAmount = nodes.size();
        return new double[nodesAmount][nodesAmount];
    }

    static double[][] agregateH(double[][] globalH, List<Element> elements) {
        for (Element element : elements) {
            double[][] elementLocalH = element.getLocalH();
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    globalH[element.getID(i)-1][element.getID(j)-1] += elementLocalH[i][j];
                }
            }
        }
        return globalH;
    }

    static double[][] agregateC(double[][] globalC, List<Element> elements) {
        for (Element element : elements) {
            double[][] elementLocalC = element.getLocalC();
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    globalC[element.getID(i)-1][element.getID(j)-1] += elementLocalC[i][j];
                }
            }
        }
        return globalC;
    }

    static void displayGlobalMatrix(double[][] globalH, List<Node> nodes) {
        int size = nodes.size();
        for (int i =0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                System.out.printf("|%-6.2f|", globalH[i][j]);
            }
            System.out.println();
        }
    }

    // ----------------------------------------- Vector P methods

    static double[] createGlobalVectorP(List<Node> nodes) {
        int nodesAmount = nodes.size();
        return new double[nodesAmount];
    }

    static double[] agregateVectorP(List<double[]> vectorsP, double[] globalVectorP, List<Element> elements) {
        int vectorCounter = 0;
        for (Element element : elements) {
            element.setVectorP(vectorsP.get(vectorCounter));
            vectorCounter++;
            double[] elementLocalVectorP = element.getVectorP();
            for (int i = 0; i < 4; i++) {
                    globalVectorP[element.getID(i)-1] += elementLocalVectorP[i];
            }
        }
        return globalVectorP;
    }

    static void displayGlobalVectorP(double[] globalVectorP, List<Node> nodes) {
        int size = nodes.size();
        for (int i =0; i < size; i++) {
            System.out.printf("|%-6.2f|", globalVectorP[i]);
        }
    }

}
