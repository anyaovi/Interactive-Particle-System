/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package predictionmodel;

import Jama.Matrix;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

/**
 *
 * @author Ana
 */
public class PredictionModel {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        ExtractSamples eS = new ExtractSamples();
        int nrSeg = 50;
        /*for running the static scenario*/
        String path = "/Users/ana.simionovici/Desktop/GOOD_EXPERIMENTS/newData/";
        String hourFirst10 = "h10";
        String hourSecond11 = "h11";
        String testPath10 = path+"h10_test.txt";
         
        String hourFirst12 = "h12";
        String hourSecond13 = "h13";      
        String testPath12 = path +"h12_test.txt";
        
        String hourFirst20 = "h20";
        String hourSecond21 = "h21";      
        String testPath20 = path +"h20_test.txt";
      
        String shuffled_hourFirst10 = "shuffled_h10";
        String  shuffled_hourSecond11 = "shuffled_h11";
        String  shuffled_testPath10 = path+"shuffled_h10_test.txt";
         
        String  shuffled_hourFirst12 = "shuffled_h12";
        String  shuffled_hourSecond13 = "shuffled_h13";      
        String  shuffled_testPath12 = path +"shuffled_h12_test.txt";
        
        String  shuffled_hourFirst20 = "shuffled_h20";
        String  shuffled_hourSecond21 = "shuffled_h21";      
        String  shuffled_testPath20 = path +"shuffled_h20_test.txt";
        /*d- dimension of the space; n - degrees of freedom*/
        int d = 2;
        int n = 5;

        int steps = 100;
        double a = 0.9999;
        int q = 1;      
          /*percentage of the particles that will survive during pertubation*/
        int percentage = 20;
        
        int iterations = 30;
        double[] z1 = eS.readTestSet(shuffled_testPath10);
        double[] z2 = eS.readTestSet(shuffled_testPath12);
        double[] z3 = eS.readTestSet(shuffled_testPath20);
        
      for(int i = 0 ; i < iterations; i++)
        {
             runAlgorithmStatic( i, eS, d, n, a, steps, percentage, q, path,  shuffled_hourFirst10,  shuffled_hourSecond11, z1);  
             runAlgorithmStatic( i, eS, d, n, a, steps, percentage, q, path,  shuffled_hourFirst12,  shuffled_hourSecond13, z2);  
             runAlgorithmStatic( i, eS, d, n, a, steps, percentage, q, path,  shuffled_hourFirst20,  shuffled_hourSecond21, z3);      
        }
    
        String pathAll = "/Users/ana.simionovici/Desktop/GOOD_EXPERIMENTS/newData/since2013";
        String hour10 = "H10";
        String hour11 = "H11";

        String hour12 = "H12";
        String hour13 = "H13";

        String hour20 = "H20";
        String hour21 = "H21";

        int windowsize = 594;

     
  /*     for(int i = 0 ; i < iterations; i++)
    {
         runAlgorithmDynamic(windowsize, i, eS, d, n, a, steps, percentage, q, pathAll, hour20, hour21);      
  }*/


    }
    /*the algorithm runs for the static scenario*/
    public static void runAlgorithmStatic(int index, ExtractSamples eS, int d, int n, double a, int steps, int percentage, int q, String path, String hourFirst, String hourSecond, double[] z) {
        String extension = "_tr.txt";
        String firstTimeFrame = path+hourFirst+extension;
        String secondTimeFrame = path+hourSecond+extension;
        double[][] samples = eS.readSamples(firstTimeFrame, secondTimeFrame);
        Matrix M = new Matrix(samples);
        double[][] samplesT = M.transpose().getArrayCopy();
        Particle init = new Particle();
        double[][] scaleMatrix = init.scaleMatrix(d);
       
    //    System.out.println("blabla samples");
      //  init.readMatrix(samples);
        
        double[] mu = init.generateMu(d);
        ArrayList<Particle> particles = generateParticles(d, mu, scaleMatrix, n, 1000, samplesT);
            //readParticlesParameters(particles);           

        //perturbation of the particles: matrix x, sigma, vector mu, recalculate L for number of steps = steps
        do {
            pert(particles, a, d, n, scaleMatrix, samplesT, percentage, steps);
            selection(particles, percentage);
            steps--;
        } while (steps > 0);

            double[] result = averageForExtractZ(particles, z, q);
            String result_Average_Path = "/Users/ana.simionovici/Desktop/GOOD_EXPERIMENTS/newData/output/resultAv" + hourSecond + "_" + index + ".txt";
            eS.writeResults(result, result_Average_Path);
            
            Particle P = getMaxParticle(particles);
            System.out.println("Sigma particle maxim likelihood" + P.L);
            double[] result2 = extractZ(P, 1, z);
            String result_Max_Path = "/Users/ana.simionovici/Desktop/GOOD_EXPERIMENTS/newData/output/resultMax"+hourSecond+"_" + index + ".txt";
            eS.writeResults(result2, result_Max_Path);
    
        
    }
   
    /*the algortihm runs for the dynamic scenario*/
    public static void runAlgorithmDynamic(int windowsize, int index, ExtractSamples eS, int d, int n, double a, int steps, int percentage, int q, String path, String hourFirst, String hourSecond) {
        String extension = ".txt";
        String firstTimeFrame = path + hourFirst + extension;
        String secondTimeFrame = path + hourSecond + extension;
        
        int poz_init = 0;
        int poz_fin = poz_init + windowsize;
       
        double[] z = new double[1];
        double[][] all = eS.readSamples(firstTimeFrame, secondTimeFrame);
        Particle init = new Particle();
        double[][] scaleMatrix = init.scaleMatrix(d);
      
        System.out.println("length of prediction    "+ (all.length - windowsize));
        double[] predictionMax = new double[all.length - windowsize];
        double[] predictionAv = new double[all.length - windowsize]; 
      
        do {
            int stepsPert = steps;
     //       System.out.println("windowsize  " + windowsize);
     //       System.out.println("nr columns      " + all[0].length);
            double[][] samples = new double[windowsize][all[0].length];
            int aux = 0;
            for (int i = poz_init; i < poz_fin; i++) {
                for (int j = 0; j < all[0].length; j++) {
                    samples[aux][j] = all[i][j];
                }
                aux++;
            }
            z[0] = all[poz_fin][0];
      //      System.out.println("z  pentru testare      " + z[0]);

            
            Matrix M = new Matrix(samples);
            double[][] samplesT = M.transpose().getArrayCopy();
            double[] mu = init.generateMu(d);
            ArrayList<Particle> particles = generateParticles(d, mu, scaleMatrix, n, 1000, samplesT);
            
            do {
                pert(particles, a, d, n, scaleMatrix, samplesT, percentage, steps);
                selection(particles, percentage);
                stepsPert--;
            } while (stepsPert > 0);
            System.out.println("poz for pred    "+poz_init);
            
            predictionAv[poz_init] =  averageForExtractZ(particles, z, q)[0];
            //       System.out.println("prediction      "+predictionAv[poz_init]);
            Particle P = getMaxParticle(particles);
            predictionMax[poz_init] = extractZ(P, q, z)[0];
            System.out.println("max particle    prediction  " + predictionMax[poz_init]);
            
            poz_init++;
            poz_fin++;
                    
            
        } while (poz_fin < all.length);
        String result_Average_Path = "/Users/ana.simionovici/Desktop/newExperiments/newData/output/resultAv_Window" + hourSecond + "_" + index + ".txt";
        eS.writeResults(predictionMax, result_Average_Path);
            
       String result_Max_Path = "/Users/ana.simionovici/Desktop/newExperiments/newData/output/resultMax_Window"+hourSecond+"_"+index+".txt";
       eS.writeResults(predictionAv, result_Max_Path);
    }
    /*read the parameters of the particles*/
    public static void readParticlesParameters(ArrayList<Particle> particles) {
        for (int i = 0; i < particles.size(); i++) {
            Particle G = particles.get(i);
            System.out.println("sigma   " + i);
            G.readMatrix(G.getSigma());
            System.out.println("mu  " + i);
            G.readVector(G.getMu());
            System.out.println("L   " + G.L+"   poz     "+i);
        }
    }
    /*extract particle with maximum lkelihood from a list of particles*/
    public static Particle getMaxParticle(ArrayList<Particle> particles) {
        ArrayList<Double> L = new ArrayList();
        for (int i = 0; i < particles.size(); i++) {
            L.add(particles.get(i).getL());
        }
        int pos = L.indexOf(Collections.max(L));
  //      System.out.println(" max particle L pos " + pos );
        return particles.get(pos);
    }
    /*extract particle with minimum lkelihood from a list of particles*/
    public static Particle getMinParticle(ArrayList<Particle> particles) {
        ArrayList<Double> L = new ArrayList();
        for (int i = 0; i < particles.size(); i++) {
            L.add(particles.get(i).getL());
        }
        int pos = L.indexOf(Collections.min(L));
        System.out.println("  min L  pos " + pos);
        return particles.get(pos);
    }
    /*generation of Particles; set up of encoded parameters*/
    public static ArrayList<Particle> generateParticles(int d, double[] mu, double[][] scaleMatrix, int n, int size, double[][] samples) {
        ArrayList<Particle> particles = new ArrayList();
        for (int i = 0; i < size; i++) {
            Particle P = new Particle();
   
            double[] muDiff = {0.0, 0.0};
            double[][] X = P.generateParticleVector(d, scaleMatrix, muDiff, n);
            double[][] sigma = P.calculateSigma(X, n);
            double L = P.calculL(samples, sigma, mu);
            P.setL(L);
            P.setSigma(sigma);
            P.setX(X);
            P.setMu(mu);
            particles.add(P);
        }
        return particles;
    }
    /*perturbation of the particles*/
    public static void pert(ArrayList<Particle> particles, double a, int d, int n, double[][] scaleMatrix, double[][] samples, int percentage, int aux) {

        for (int i = 0; i < particles.size(); i++) {
            Particle p = particles.get(i);
            double[][] x = p.getX();
            double L = p.getL();
            double[] mu = p.getMu();
            double[][] newX = x;
            double[] muDiff = {0.0, 0.0};
            double[][] y = p.generateParticleVector(d, scaleMatrix, muDiff, n);            
            newX = p.perturbX(newX, y, a, d, n);
            double[] muNew = p.perturbMu(mu, aux);
            double[][] newsigma = p.calculateSigma(newX, n);
     //       System.out.println("Particle    "+ i+"   L before    " +L);
            double Lnew = p.calculL(samples, newsigma, muNew);
      //      System.out.println("L perturbat     "+Lnew );
//tranzitie
            if (Lnew > L) {
                p.setX(newX);
                p.setSigma(newsigma);
                p.setL(Lnew);
                p.setMu(muNew);
            }
          //   System.out.println("L new is    "+ p.getL()+"\n");
        }           
  //      readParticlesParameters(particles);
        
    
    }
    /*selection of a percentage of particles based on the likelihood*/
    public static void selection(ArrayList<Particle> particles, int percentage)    {
        int size = particles.size();
        int positionBound = size*percentage/100;          
 //       System.out.println("position  bound  " +positionBound);
        
        ArrayList<Double> L = new ArrayList();
        for (int i = 0; i < particles.size(); i++) {
            L.add(particles.get(i).getL());
        }        
        Collections.sort(L, Collections.reverseOrder());
             /*   for(int i = 0 ; i < L.size();i++)
        {
            System.out.println(i+" "+L.get(i));
        }   
            System.out.println("");
       */ 
        double Lbound = L.get(positionBound);
   //     System.out.println("L to compare with       "+Lbound);        
        ArrayList<Particle> selected = new ArrayList(positionBound);     
        
        for(int i = 0 ; i < particles.size();i++)
        {
            Particle P = particles.get(i);
            if(P.getL() >= Lbound)
            {
                selected.add(P);
           //     System.out.println("heheh       L   "+P.getL());
            }
        }
        
      
      for(int i = 0 ; i <size;i++)
        {
            Particle P = particles.get(i);
            if(P.getL() < Lbound)
            {
                Random rnd = new Random();
                int pos = rnd.nextInt(positionBound);
                Particle M = selected.get(pos);
                P.setL(M.getL());
                P.setMu(M.getMu());
                P.setSigma(M.getSigma());
                P.setX(M.getX());             
             }
        }
             
       // readParticlesParameters(particles);
     }
    /*prediction method using the particle with the maximum likelihood from the final population*/
    public static double[] extractZ(Particle P, int q, double[] z) {

        double[] mu = P.getMu();
        double[][] sigma = P.getSigma();
        int n = P.getX()[0].length;
        int d = P.getX().length;

   //     System.out.println("q is    " + q);
        double dim = d - q;
   //     System.out.println("d-q     is " + dim);

        double[] mualpha = new double[q];
        double[] mubeta = new double[d - q];
        double[][] sigmaAlpha1 = new double[q][q];
        double[][] sigmaAlpha2 = new double[q][d - q];
        double[][] sigmaBeta1 = new double[d - q][q];
        double[][] sigmaBeta2 = new double[d - q][d - q];

        double[] muC = new double[q];
        double[][] sigmaC = new double[q][q];
        System.arraycopy(mu, 0, mualpha, 0, q);
        System.arraycopy(mu, q, mubeta, 0, d - q);
   //     System.out.println("mu alpha    " + mualpha[0]);
   //     System.out.println("mu beta     " + mubeta[0]);
   //     System.out.println("SIGMA");
   //     P.readMatrix(sigma);

        for (int i = 0; i < q; i++) {
            System.arraycopy(sigma[i], 0, sigmaAlpha1[i], 0, q);
        }

        for (int i = 0; i < q; i++) {
            for (int j = 0; j < d - q; j++) {
                sigmaAlpha2[i][j] = sigma[i][j + q];
            }
        }

        for (int i = 0; i < d - q; i++) {
            System.arraycopy(sigma[i + q], 0, sigmaBeta1[i], 0, q);
            for (int j = 0; j < d - q; j++) {
                sigmaBeta2[i][j] = sigma[i + q][j + q];
            }
        }

        /*System.out.println("MU  ");
         readVector(mu);
         System.out.println("mualpha");
         readVector(mualpha);
         System.out.println("mubeta");
         readVector(mubeta);
         System.out.println("sigma");
         P.afisareM(sigma);
         System.out.println("sigma alfa 1");
         P.afisareM(sigmaAlpha1);
         System.out.println("sigma alfa 2");
         P.afisareM(sigmaAlpha2);
         System.out.println("sigma beta 1");
         P.afisareM(sigmaBeta1);
         System.out.println("sigma beta 2");
         P.afisareM(sigmaBeta2);
         */
        Matrix M = new Matrix(sigmaAlpha1);
        double[][] sigmaAlpha1Inverse = M.inverse().getArrayCopy();
        /*   System.out.println("Inverse for alpha1");
         P.afisareM(sigmaAlpha1Inverse);
       
         System.out.println("z           "+z[0]);*/

        double[] substraction = new double[z.length];

        for (int i = 0; i < z.length; i++) {
            substraction[i] = z[i] - mualpha[0];
        }
       
        double[][] productM = P.multiplyMM(sigmaBeta1, sigmaAlpha1Inverse);        
        double[] productMV = P.multiplyMV(productM, substraction);      
        double[] zbeta = new double[productMV.length];
        for (int i = 0; i < productMV.length; i++) {
            zbeta[i] = mubeta[0] + productMV[i];
        }

        /*        System.out.println("zbeta   length" + zbeta.length);
         readVector(zbeta);
         */
        /*  old version X2 X1
         Matrix M = new Matrix(sigmaBeta2);
         double[][] sigmaBeta2Inverse = M.inverse().getArrayCopy();
         System.out.println("Inverse for beta2");
         P.afisareM(sigmaBeta2Inverse);

         double[] substraction = new double[z.length];

         for (int i = 0; i < z.length; i++) {
         substraction[i] = z[i] - mubeta[0];
         }
         System.out.println("substraction");
         readVector(substraction);

         double[][] productM = P.multiply(sigmaAlpha2, sigmaBeta2Inverse);
         System.out.println("productM  alpha2 cu beta2 inverse");
         P.afisareM(productM);

         double[] productMV = P.multiplyMV(productM, substraction);
         System.out.println("productM cu substraction");
         readVector(productMV);
    
       
         double[] zbeta = new double[productMV.length];
         for (int i = 0; i < productMV.length; i++) {
         zbeta[i] = mualpha[0] + productMV[i];
         }

         System.out.println("zbeta   length" + zbeta.length);
         readVector(zbeta);

      
         //   System.out.println("Max likelihood");
         //  P.afisareM(result);
         // System.out.println("Max Likelihood ------------ just result");

         //   readVector(zbeta);
         */
        return zbeta;

    }
    /*prediction method using the average sum of the likelihoods of the particles from the final population*/
    public static double[] averageForExtractZ(ArrayList<Particle> particles, double[] z, int q) {
        double sumL = 0.0;
        double[] zbeta = new double[z.length];
        double[] L = new double[particles.size()];
        for (int i = 0; i < particles.size(); i++) {
            sumL += particles.get(i).getL();
            L[i] = particles.get(i).getL();
        }
        //
        double[][] predictions = new double[z.length][particles.size()];
        for (int i = 0; i < particles.size(); i++) {
            Particle P = particles.get(i);
            double[] pred = extractZ(P, q, z);
            for (int j = 0; j < pred.length; j++) {
                predictions[j][i] = pred[j];
            }
       }
//calculeaza average
      //   System.out.println("L is");
    //     particles.get(0).readVector(L);
   //      System.out.println("Prediction");
     //    particles.get(0).readMatrix(predictions);

       for(int i = 0; i < zbeta.length;i++)
       {
           double val = 0.0; 
           for(int j = 0; j < particles.size();j++)
           {
               double l = L[j];
               val += l * predictions[i][j];
           }
           zbeta[i] = val/sumL;
       }
         
        //System.out.println("average extraction");
        //  particles.get(0).afisareM(result);
        System.out.println("average extraction ----------------- just result");
        particles.get(0).readVector(zbeta);
        return zbeta;

    }
}
