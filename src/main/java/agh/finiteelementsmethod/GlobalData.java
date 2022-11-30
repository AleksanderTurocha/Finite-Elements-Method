package agh.finiteelementsmethod;

import lombok.Data;

@Data
public class GlobalData {
    private double simulationTime;
    private double simulationStepTime;
    private double conductivity;
    private double alfa;
    private double tot;
    private double initialTemp;
    private double density;
    private double specificHeat;

    public GlobalData(double simulationTime, double simulationStepTime, double conductivity, double alfa, double tot, double initialTemp, double density, double specificHeat) {
        this.simulationTime = simulationTime;
        this.simulationStepTime = simulationStepTime;
        this.conductivity = conductivity;
        this.alfa = alfa;
        this.tot = tot;
        this.initialTemp = initialTemp;
        this.density = density;
        this.specificHeat = specificHeat;
    }
}
