package agh.finiteelementsmethod;

import java.util.List;

// To change points you have to change it here in mesapp and also if you want to work with certain agregated matrix you have to change points in BorderConditions -> calculateBorderConditions
public class MesApp {
    public static void main(String[] args) {
        //Basics of FEM
        MeshReader meshReader = new MeshReader();
        List<String> fileLines = meshReader.readData();
        GlobalData globalData = meshReader.createGlobalData(fileLines);
        List<Node> nodes = meshReader.createNodes(fileLines);

        List<Element> elements = meshReader.createElements((fileLines));
        meshReader.assignEdgeConditions(nodes, meshReader.createEdgeConditions(fileLines));
        Grid grid = meshReader.createGrid(fileLines, nodes, elements);

//        System.out.println(fileLines);
//        System.out.println(globalData);
//        System.out.println(nodes);
//        System.out.println(elements);
//        System.out.println(globalData);

        //Integration
//        Integration integration = new Integration();
//        System.out.println(integration.integrate1D(2));
//        System.out.println(integration.integrate1D(3));
//        System.out.println(integration.integrate2D(2));
//        System.out.println(integration.integrate2D(3));

        //Matrix H
        HMatrix hMatrix = new HMatrix(globalData, grid, nodes, elements);
        double[][] H2P = hMatrix.calculateMatrixH(2, 0);

//        // 2P
//        System.out.println("2 points:");
//        double[][] H2P = hMatrix.calculatingMatrixH(2);
//        hMatrix.displayH(H2P);
//
//        // 3P
//        System.out.println("3 points:");
//        double[][] H3P = hMatrix.calculatingMatrixH(3);
//        hMatrix.displayH(H3P);
//
        // 4P
//        System.out.println("4 points:");
//        double[][] H4P = hMatrix.calculatingMatrixH(4);
//        hMatrix.displayH(H4P);

        // Matrix Hbc
        BorderConditions borderConditions = new BorderConditions(nodes, elements, globalData, new UniversalElement(), grid, hMatrix);

        // Agregation with border conditions
        System.out.println("Here starts global H with border conditions:");
        double[][] borderConditionsMatrix = borderConditions.calculateBorderConditions();
        Agregation.displayGlobalH(borderConditionsMatrix, nodes);
    }
}