package agh.finiteelementsmethod;

public class Point {
    private double ksi;
    private double eta;

    public double getKsi() {
        return ksi;
    }

    public void setKsi(double ksi) {
        this.ksi = ksi;
    }

    public double getEta() {
        return eta;
    }

    public void setEta(double eta) {
        this.eta = eta;
    }

    public Point(double ksi, double eta) {
        this.ksi = ksi;
        this.eta = eta;
    }
}
