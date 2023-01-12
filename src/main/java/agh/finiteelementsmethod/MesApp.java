package agh.finiteelementsmethod;

import java.util.List;

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
//        double[][] H2P = hMatrix.calculateMatrixH(2, 0);

//        // 2P
//        System.out.println("2 points:");
//        double[][] H2P = hMatrix.calculateMatrixH(2, 0);
//        hMatrix.displayH(H2P);
//
//        // 3P
//        System.out.println("3 points:");
//        double[][] H3P = hMatrix.calculateMatrixH(3, 0);
//        hMatrix.displayH(H3P);
//
//        // 4P
//        System.out.println("4 points:");
//        double[][] H4P = hMatrix.calculateMatrixH(4, 0);
//        hMatrix.displayH(H4P);

        // Matrix Hbc
        BorderConditions borderConditions = new BorderConditions(nodes, elements, globalData, new UniversalElement(), grid, hMatrix);

        // Agregation matrix H with border conditions
        System.out.println("Here starts global H with border conditions:");
        double[][] borderConditionsMatrixH = borderConditions.calculateBorderConditions(2);
        Agregation.displayGlobalMatrix(borderConditionsMatrixH, nodes);

        // Vector P with agregation
        System.out.println("Here starts vector P with agregation:");
        double[] globalVectorP = borderConditions.calculateVectorP(2);
        Agregation.displayGlobalVectorP(globalVectorP, nodes);

        // Matrix C
        System.out.println();
        System.out.println("Here starts global matrix C:");
        CMatrix cMatrix = new CMatrix(globalData, grid, nodes, elements);
        double[][] finalMatrixC = borderConditions.calculateMatrixC(2);
        Agregation.displayGlobalMatrix(finalMatrixC, nodes);

        // Time
        Time.applyTime(borderConditionsMatrixH, finalMatrixC, globalVectorP, globalData, nodes);

    }
}