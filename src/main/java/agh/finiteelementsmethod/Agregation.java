package agh.finiteelementsmethod;

import java.util.List;

public class Agregation {
    static double[][] createGlobalH(List<Node> nodes) {
        int nodesAmount = nodes.size();
        return new double[nodesAmount][nodesAmount];
    }

    static double[][] agregate(double[][] globalH, List<Element> elements) {
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

    static void displayGlobalH(double[][] globalH, List<Node> nodes) {
        int size = nodes.size();
        for (int i =0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                System.out.printf("|%-6.2f|", globalH[i][j]);
            }
            System.out.println();
        }
    }
}
