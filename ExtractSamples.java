/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package predictionmodel;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author Ana
 */
public class ExtractSamples {

    public ExtractSamples() {
    }

    public int[][] count(ArrayList<Integer> allfirst, ArrayList<Integer> allsecond, int nrSeg, boolean cut) {
        ArrayList<Days> allDays = new ArrayList();
        double av1 = 0.0;
        double av2 = 0.0;
        for (int i = 0; i < allfirst.size(); i++) {
            av1 += allfirst.get(i);
            av2 += allsecond.get(i);
        }
        av1 = av1 / allfirst.size();
        av2 = av2 / allfirst.size();

        int maxX1 = Collections.max(allfirst);
        int maxX2 = Collections.max(allsecond);

        //adding this condition we avoid adding the days with a smaller number of calls, for hour 10 and 11, then average
        if (cut == true) {
            for (int i = 0; i < allfirst.size(); i++) {
                if (allfirst.get(i) >= av1 && allsecond.get(i) >= av2) {
                    Days d = new Days();
                    d.setX1(allfirst.get(i));
                    d.setX2(allsecond.get(i));
                    allDays.add(d);
                }
            }
        }
        else{

//this adds all days from the database no matter what the number of calls per hour 10 and 11 is
        for (int i = 0; i < allfirst.size(); i++) {
                Days d = new Days();          
                d.setX1(allfirst.get(i));
                d.setX2(allsecond.get(i));
                allDays.add(d);

            }
        }
        ArrayList<Integer> axisX1 = new ArrayList();
        ArrayList<Integer> axisX2 = new ArrayList();

        System.out.println("nr days " + allDays.size());


        int[][] countDays = new int[nrSeg][nrSeg];
        int sizeIntX1 = Math.round(maxX1 / nrSeg);
        //      System.out.println("size interval x1    " + sizeIntX1);
        int sizeIntX2 = Math.round(maxX2 / nrSeg);
        //     System.out.println("size interval x2    " + sizeIntX2);

        for (int i = 0; i < allDays.size(); i++) {
            //       System.out.println("");
            Days d = allDays.get(i);
            //       System.out.println("data    " + i);



            for (int k = 0; k < nrSeg; k++) {

                if ((d.getX1() > (k * sizeIntX1)) && (d.getX1() < ((k + 1) * sizeIntX1)) || (d.getX1() == ((k + 1) * sizeIntX1))) {
//                    System.out.println("x1     "
//                            + d.getX1() + "  step  " + k + " k * sizeIntX1    " + k * sizeIntX1
//                            + "  ( k + 1) * sizeIntX1    " + (k + 1) * sizeIntX1);


                    for (int j = 0; j < nrSeg; j++) {

                        if ((d.getX2() > (j * sizeIntX2)) && (d.getX2() < ((j + 1) * sizeIntX2)) || (d.getX2() == ((j + 1) * sizeIntX2))) {

//                            System.out.println("x2     "
//                                    + d.getX2() + "  step  " + j + " j * sizeIntX2    "
//                                    + j * sizeIntX2 + "  ( j + 1) * sizeIntX2    " + (j + 1) * sizeIntX2);
                            countDays[k][j]++;
                            //  System.out.println("   x1  "+k+"  x2 " +j+"    count   " +countDays[k][j]);
                            break;

                        }
                        if (d.getX2() > (nrSeg * sizeIntX2)) {
                            countDays[k][nrSeg - 1]++;
//                            System.out.println("x2     "
//                                    + d.getX2() + "  step  " + j + " j * sizeIntX2    "
//                                    + (nrSeg - 1) + "  ( j + 1) * sizeIntX2    " + nrSeg * sizeIntX2);
                            break;
                        }
                    }

                } else if (d.getX1() > (nrSeg * sizeIntX1)) {
                    for (int j = 0; j < nrSeg; j++) {

                        if ((d.getX2() > (j * sizeIntX2)) && (d.getX2() < ((j + 1) * sizeIntX2)) || (d.getX2() == ((j + 1) * sizeIntX2))) {
//
//                            System.out.println("x2     "
//                                    + d.getX2() + "  step  " + j + " j * sizeIntX2    "
//                                    + j * sizeIntX2 + "  ( j + 1) * sizeIntX2    " + (j + 1) * sizeIntX2);
                            countDays[nrSeg - 1][j]++;
                            // System.out.println("   x1  " + k + "  x2 " + j + "    count   " + countDays[nrSeg - 1][j]);
                            break;

                        }
                        if (d.getX2() > (nrSeg * sizeIntX2)) {

                            countDays[nrSeg - 1][nrSeg - 1]++;
                            break;
                        }
                    }

                    break;
                }
            }
        }
        for (int i = 0; i < nrSeg; i++) {
            axisX1.add(i * sizeIntX1);
            axisX2.add(i * sizeIntX2);
        }
        for (int i = 0; i < axisX1.size(); i++) {
            System.out.print(axisX1.get(i) + ", ");


        }
        System.out.println();
        for (int j = 0; j < axisX2.size(); j++) {
            System.out.print(+axisX2.get(j) + ", ");
        }
        System.out.println();
        return countDays;
    }

    public int sumElem(int[][] matrix) {
        int sum = 0;
        for (int i = 0; i < matrix[0].length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                sum += matrix[i][j];
            }
        }
        return sum;
    }   
   
    public ArrayList<Integer> readFile(String path) {
        ArrayList<Integer> values = new ArrayList();

        //Name of the file
      //  String fileName = path;
        try {

            //Create object of FileReader
            FileReader inputFile = new FileReader(path);

            //Instantiate the BufferedReader Class
            BufferedReader bufferReader = new BufferedReader(inputFile);

            //Variable to hold the one line data
            String line;

            // Read file line by line and print on the console
            while ((line = bufferReader.readLine()) != null) {
                //   System.out.println(line);
                values.add(Integer.parseInt(line));
                // System.out.println("");
            }
            //  System.out.println("lenght  " + values.size() + "    first value " + values.get(2));
            //Close the buffer reader
            bufferReader.close();
        } catch (Exception e) {
            System.out.println("Error while reading file line by line:"
                    + e.getMessage());
        }
        return values;
    }

    public double[] readTestSet(String path) {
       ArrayList<Integer> testset = readFile(path);

        double[] z = new double[testset.size()];
        for(int i = 0; i< testset.size(); i++)
        {
            z[i] = testset.get(i);
        }
            return z;
    }
    
    public void writeResults(double[] z, String path)
    {
        try {            //Create object of FileWriter
           FileWriter fw = new FileWriter(new File(path));           
           for(int i = 0 ; i < z.length; i++)
           {
               fw.write(String.valueOf(z[i]));
               fw.write(System.lineSeparator());
           }                      
            fw.close();
        } catch (Exception e) {
            System.out.println("Error while writing file"
                    + e.getMessage());
        }
    }
    
    public double[][] calcSamples(int[][] count) {
        int sumElem = sumElem(count);
        double[][] result = new double[sumElem][2];
        int pos = 0;
        for (int i = 0; i < count.length; i++) {
            for (int j = 0; j < count[0].length; j++) {
                int val = count[i][j];

                if (val > 0) {
                    for (int k = pos; k < pos + val; k++) {
                        result[k][0] = i;
                        result[k][1] = j;
                    }
                    pos = pos + val;
                }
            }
        }
        return result;
    }

    public double[] calcMean(double[][] samples) {
        double[] mean = new double[samples[0].length];
        for (int i = 0; i < samples[0].length; i++) {
            double sum = 0.0;
            for (int j = 0; j < samples.length; j++) {
                sum += samples[j][i];
            }
            mean[i] = sum / samples.length;
            System.out.println("sum     " + sum + "  mean    " + " i   " + i + "   " + mean[i]);
        }

        return mean;
    }

    public double[][] addScalarMatrix(double[][] matrix, double a) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] = matrix[i][j] + a;
            }
        }
        return matrix;
    }

    public double[][] readSamples(String first, String second) {
        ArrayList<Integer> allfirst = readFile(first);
        System.out.println("allfirst size   " + allfirst.size());
        ArrayList<Integer> allsecond = readFile(second);
      //   int[][] countDays = count(allfirst, allsecond, nrSeg,false);
      //   afisareM(countDays);
      //   double[][] samples = calcSamples(countDays);
        double[][] samples = new double[allfirst.size()][2];
        
       for(int i = 0 ; i < allfirst.size();i++)
               {
                    samples[i][0] = allfirst.get(i);
                    samples[i][1] = allsecond.get(i);
                }
                return samples;
            }
      
    public int[][] createCountMatrix(double[][] result) {

        ArrayList<Double> allValues = new ArrayList();
        for (int i = 0; i < result.length; i++) {

            allValues.add(result[i][0]);

        }
        int nrSeg = Collections.max(allValues).intValue() + 1;
        System.out.println("nrSeg   " + nrSeg);
        int[][] countDays = new int[nrSeg][nrSeg];


        for (int i = 0; i < result.length; i++) {
            Double x1 = result[i][0];
            //   System.out.println("x1  " + x1);
            Double x2 = result[i][1];
            //     System.out.println("x2  " + x2);
            for (int j = 0; j < nrSeg; j++) {
                if (j == x1.intValue()) {
                    for (int k = 0; k < nrSeg; k++) {
                        if (k == x2.intValue()) {
                            countDays[j][k] = countDays[j][k] + 1;
                        }
                    }
                }
            }
        }
        return countDays;

    }
    
   
}
