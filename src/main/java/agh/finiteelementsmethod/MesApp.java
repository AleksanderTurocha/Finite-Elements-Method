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

        System.out.println(fileLines);
        System.out.println(globalData);
        System.out.println(nodes);
        System.out.println(elements);
        System.out.println(globalData);

        //Integration
//        Integration integration = new Integration();
//        System.out.println(integration.integrate1D(2));
//        System.out.println(integration.integrate1D(3));
//        System.out.println(integration.integrate2D(2));
//        System.out.println(integration.integrate2D(3));

        //Universal elements

    }
}