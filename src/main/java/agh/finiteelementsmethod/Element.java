package agh.finiteelementsmethod;

import java.util.Arrays;

public class Element {
    private final int elementID;
    private final int[] ID = new int[4];

    private double[][] localH;

    private double[] vectorP;

    private double[][] localC;

    public Element(int elementID, int n1, int n2, int n3, int n4) {
        this.elementID = elementID;
        this.ID[0] = n1;
        this.ID[1] = n2;
        this.ID[2] = n3;
        this.ID[3] = n4;
    }

    public double[][] getLocalC() {
        return localC;
    }

    public void setLocalC(double[][] localC) {
        this.localC = localC;
    }

    public double[] getVectorP() {
        return vectorP;
    }

    public void setVectorP(double[] vectorP) {
        this.vectorP = vectorP;
    }

    public double[][] getLocalH() {
        return localH;
    }

    public void setLocalH(double[][] localH) {
        this.localH = localH;
    }

    public int getID(int index) {
        return ID[index];
    }

    @Override
    public String toString() {
        return "Element{" +
                "ID=" + Arrays.toString(ID) +
                '}';
    }
}
