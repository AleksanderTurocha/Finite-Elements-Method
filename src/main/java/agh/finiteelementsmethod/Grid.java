package agh.finiteelementsmethod;

import lombok.Data;

import java.util.ArrayList;
import java.util.List;

@Data
public class Grid {
    private int nodesNumber;
    private int elementsNumber;
    private List<Node> nodes = new ArrayList<>();
    private List<Element> elements = new ArrayList<>();

    public Grid(int nodesNumber, int elementsNumber, List<Node> nodes, List<Element> elements) {
        this.nodesNumber = nodesNumber;
        this.elementsNumber = elementsNumber;
        this.nodes = nodes;
        this.elements = elements;
    }
}
