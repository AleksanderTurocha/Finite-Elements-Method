package agh.finiteelementsmethod;

import java.util.Arrays;

public class Element {
    private int[] ID = new int[4];

    public Element(int n1, int n2, int n3, int n4) {
        this.ID[0] = n1;
        this.ID[1] = n2;
        this.ID[2] = n3;
        this.ID[3] = n4;
    }

    @Override
    public String toString() {
        return "Element{" +
                "ID=" + Arrays.toString(ID) +
                '}';
    }
}
