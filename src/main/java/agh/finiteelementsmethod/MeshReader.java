package agh.finiteelementsmethod;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

public class MeshReader {

    List<String> readData() {
        final String fileName = "Mesh.txt";
        try {
            BufferedReader reader = new BufferedReader(new FileReader(fileName));
            return reader.lines().toList();
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    GlobalData createGlobalData(List<String> fileLines) {
        int iterator = 0;
        List<Double> globalDataNumbers = new ArrayList<>();
        while (!fileLines.get(iterator).equals("*Node")) {
            globalDataNumbers.add(Double.parseDouble(fileLines.get(iterator).split(" ")[1]));
            iterator++;
        }
        return new GlobalData(
            globalDataNumbers.get(0),
            globalDataNumbers.get(1),
            globalDataNumbers.get(2),
            globalDataNumbers.get(3),
            globalDataNumbers.get(4),
            globalDataNumbers.get(5),
            globalDataNumbers.get(6),
            globalDataNumbers.get(7)
        );
    }

    List<Node> createNodes(List<String> fileLines) {
        int iterator = 0;
        while (!fileLines.get(iterator).equals("*Node"))
            iterator++;
        iterator++;
        List<Node> nodes = new ArrayList<>();
        while (!fileLines.get(iterator).equals("*Element, type=DC2D4")) {
            nodes.add(new Node(
                    Double.parseDouble(fileLines.get(iterator).trim().split(", ")[1]),
                    Double.parseDouble(fileLines.get(iterator).trim().split(", ")[2]),
                    0,
                    0,
                    Integer.parseInt(fileLines.get(iterator).trim().split(", ")[0])
            ));
            iterator++;
        }
        assignEdgeConditions(nodes, createEdgeConditions(fileLines));
        return nodes;
    }

    int[] createEdgeConditions(List<String> fileLines) {
        int iterator = 0;
        while (!fileLines.get(iterator).equals("*BC"))
            iterator++;
        iterator++;
        int length = fileLines.get(iterator).trim().split(", ").length;
        int[] conditions = new int[length];
        for (int i = 0; i < length; i++) {
            conditions[i] = Integer.parseInt(fileLines.get(iterator).trim().split(", ")[i]);
        }
        return conditions;
    }

    void assignEdgeConditions(List<Node> nodes, int[] edgeConditionsArray) {
        int size = nodes.size();
        for (int i = 0; i <= size; i++) {
            for (int condition : edgeConditionsArray) {
                if (i + 1 == condition) {
                    nodes.get(i).setBC(1);
                }
            }
        }
    }

    List<Element> createElements(List<String> fileLines) {
        int iterator = 0;
        List<Element> elements = new ArrayList<>();
        while (!fileLines.get(iterator).equals("*Element, type=DC2D4"))
            iterator++;
        iterator++;
        while (!fileLines.get(iterator).equals("*BC")) {
            elements.add(new Element(
                    (int) Double.parseDouble(fileLines.get(iterator).trim().split(", ")[0]),
                    (int) Double.parseDouble(fileLines.get(iterator).trim().split(", ")[1]),
                    (int) Double.parseDouble(fileLines.get(iterator).trim().split(", ")[2]),
                    (int) Double.parseDouble(fileLines.get(iterator).trim().split(", ")[3]),
                    (int) Double.parseDouble(fileLines.get(iterator).trim().split(", ")[4])
            ));
            iterator++;
        }
        return elements;
    }

    Grid createGrid(List<String> fileLines, List<Node> nodeList, List<Element> elementList) {
        int iterator = 0;
        List<Integer> globalDataNumbers = new ArrayList<>();
        while (!fileLines.get(iterator).equals("*Node")) {
            globalDataNumbers.add(Integer.parseInt(fileLines.get(iterator).split(" ")[1]));
            iterator++;
        }
        return new Grid(
                globalDataNumbers.get(8),
                globalDataNumbers.get(9),
                nodeList,
                elementList
        );
    }
}
