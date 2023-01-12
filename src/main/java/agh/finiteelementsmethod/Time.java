package agh.finiteelementsmethod;

import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;

public class Time {

    static void applyTime(double[][] borderConditionsMatrixH, double[][] finalMatrixC, double[] globalVectorP, GlobalData globalData, List<Node> nodes) {
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                borderConditionsMatrixH[i][j] += finalMatrixC[i][j] / globalData.getSimulationStepTime();
            }
        }

        // Creating copy of vector p
        double[] globalVectorPCopy = Arrays.copyOf(globalVectorP, globalVectorP.length);

        // Creating copy of matrix H
        double[][] borderConditionsMatrixHCopy = new double[nodes.size()][nodes.size()];
        for(int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                borderConditionsMatrixHCopy[i][j] = borderConditionsMatrixH[i][j];
            }
        }

        // Temperature vector
        double[] temperatureVector = new double[nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            temperatureVector[i] = globalData.getInitialTemp();
        }

        for (int i = 0; i < globalData.getSimulationTime(); i += globalData.getSimulationStepTime()) {
//            System.out.println(Arrays.toString(temperatureVector));
            for (int j = 0; j < nodes.size(); j++) {
                for (int k = 0; k < nodes.size(); k++) {
                    globalVectorP[j] += (finalMatrixC[j][k] / globalData.getSimulationStepTime()) * temperatureVector[k];
                }
            }
            temperatureVector = GaussElimination.solveEquationSystem(borderConditionsMatrixH, globalVectorP);
            OptionalDouble min = Arrays.stream(temperatureVector).min();
            OptionalDouble max = Arrays.stream(temperatureVector).max();
            System.out.println("Iteration: " + i);
            double minDouble = min.getAsDouble();
            double maxDouble = max.getAsDouble();
            System.out.printf("min: %f.2\n", minDouble);
            System.out.printf("max: %f.2\n", maxDouble);

            globalVectorP = Arrays.copyOf(globalVectorPCopy, globalVectorPCopy.length);

            // Bringing the values of H back
            for(int j = 0; j < nodes.size(); j++) {
                for (int k = 0; k < nodes.size(); k++) {
                    borderConditionsMatrixH[j][k] = borderConditionsMatrixHCopy[j][k];
                }
            }
        }
    }
}
