package agh.finiteelementsmethod;

import lombok.Data;

@Data
public class Node {
    private double x;
    private double y;
    private double temperature;
    private int BC;

    public Node(double x, double y, double temperature, int BC) {
        this.x = x;
        this.y = y;
        this.temperature = temperature;
        this.BC = BC;
    }
}
