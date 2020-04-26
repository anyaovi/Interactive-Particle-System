package predictionmodel;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 *
 * @author Ana
 */
public class BoxMuller {

    public BoxMuller() {
    }
    
     public double[] boxMuller(int N) {
        double[] Z =new double[N];
        for (int i = 0; i < N; i++) {
            Random rand = new Random();
            double u1 = rand.nextDouble();
            // System.out.println("u1 is: " +u1);           
            double u2 = rand.nextDouble();
            //System.out.println("u2 is: "+u2);
            double z = Math.pow((-2 * Math.log(u1)), 0.5) * Math.cos(2 * Math.PI * u2);
           
            Z[i] = z;            //System.out.println("z is: " + z);
        }
       
        return Z;
    }
     
}
