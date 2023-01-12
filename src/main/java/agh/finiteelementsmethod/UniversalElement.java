package agh.finiteelementsmethod;

// Universal Element structure containing arrays of dNdKsi and dNdEta
public class UniversalElement {

    private double[][] arraydNdKsi;
    private double[][] arraydNdEta;

    private double[][] shapeFunctionValuesInIntegrationPoints;

    public UniversalElement(double[][] arraydNdKsi, double[][] arraydNdEta) {
        this.arraydNdKsi = arraydNdKsi;
        this.arraydNdEta = arraydNdEta;
    }

    public UniversalElement(double[][] shapeFunctionValuesInIntegrationPoints) {
        this.shapeFunctionValuesInIntegrationPoints = shapeFunctionValuesInIntegrationPoints;
    }

    public UniversalElement() {}

    public double[][] getArraydNdKsi() {
        return arraydNdKsi;
    }

    public void setArraydNdKsi(double[][] arraydNdKsi) {
        this.arraydNdKsi = arraydNdKsi;
    }

    public double[][] getArraydNdEta() {
        return arraydNdEta;
    }

    public void setArraydNdEta(double[][] arraydNdEta) {
        this.arraydNdEta = arraydNdEta;
    }

    // Creating tables derivatives eta and ksi
    double dNdKsiFunctions(int index, double eta) {
        switch (index) {
            case 0 -> {
                return DN1DKsi(eta);
            }
            case 1 -> {
                return DN2DKsi(eta);
            }
            case 2 -> {
                return DN3DKsi(eta);
            }
            case 3 -> {
                return DN4DKsi(eta);
            }
        }
        return 0;
    }
    double dNdEtaFunctions(int index, double ksi) {
        switch (index) {
            case 0 -> {
                return DN1DEta(ksi);
            }
            case 1 -> {
                return DN2DEta(ksi);
            }
            case 2 -> {
                return DN3DEta(ksi);
            }
            case 3 -> {
                return DN4DEta(ksi);
            }
        }
        return 0;
    }

    // Shape functions
    // N1
    double N1(double ksi, double eta) { return 0.25 * (1 - ksi) * (1 - eta); }

    // N2
    double N2(double ksi, double eta) {
        return 0.25 * (1 + ksi) * (1 - eta);
    }

    // N3
    double N3(double ksi, double eta) {
        return 0.25 * (1 + ksi) * (1 + eta);
    }

    // N4
    double N4(double ksi, double eta) { return 0.25 * (1 - ksi) * (1 + eta); }

    // Derivatives of shape functions
    // N1 / ksi
    double DN1DKsi(double eta) {
            return -0.25 * (1 - eta);
    }

    // N2 / ksi
    double DN2DKsi(double eta) {
            return 0.25 * (1 - eta);
    }

    // N3 / ksi
    double DN3DKsi(double eta) {
            return 0.25 * (1 + eta);
    }

    // N4 / ksi
    double DN4DKsi(double eta) {
            return -0.25 * (1 + eta);
    }

    // N1 / eta
    double DN1DEta(double ksi) {
            return -0.25 * (1 - ksi);
    }

    // N2 / eta
    double DN2DEta(double ksi) {
            return -0.25 * (1 + ksi);
    }

    // N3 / eta
    double DN3DEta(double ksi) {
            return 0.25 * (1 + ksi);
    }

    // N4 / eta
    double DN4DEta(double ksi) {
            return 0.25 * (1 - ksi);
    }
}
