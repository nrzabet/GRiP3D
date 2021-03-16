package objects;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import utils.Utils;

/**
 * class that constructs the 3D contact matrix based on the file inputted 
 * @author adumita@essex.ac.uk
 *
 */
public class InteractionMatrix implements Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private boolean[][] SimulatedMatrix;
	private double[][] CumulativeSimulatedMatrix;
	private double[][] ContactMatrix;
    private ArrayList<DNAregion> bins;
    private double[][] NormalisedMatrix;
    private double LastTimeSimulated;

	/**
	 * class constructor
	 *
	 * @param ContactMatrixLocation the path to the file containing the HiC interactions
	 * @param regionOfInterest object of type DNAregion that contains the start and end of the sequence based on which the bins are created
	 * @param binWidth the size of the bins to be created
	 * @param generator
	 * @throws Exception
	 */
	public InteractionMatrix(String ContactMatrixLocation,DNAregion regionOfInterest, int binWidth,Random generator) throws Exception {
			
			LastTimeSimulated=0;
			bins = new ArrayList<DNAregion>();
			createMatrix(regionOfInterest,binWidth);
		    addScoreToMatrix(regionOfInterest, ContactMatrixLocation); 
		    normaliseMatrix();
		    simulateMatrix(generator,0);
	}
	
	/**
	 * class constructor
	 * 
	 * @param regionOfInterest object of type DNAregion that contains the start and end of the sequence based on which the bins are created
	 * @param binWidth the size of the bins to be created
	 * @param generator
	 */
	public InteractionMatrix(DNAregion regionOfInterest,int binWidth,Random generator) {
			bins = new ArrayList<DNAregion>();
			createMatrix(regionOfInterest,binWidth);
			normaliseMatrix();
			simulateMatrix(generator,0);
	}
	
		/**
		 * creates the empty matrix using the DNA region of interest 
		 * @param regionOfInterest the subsequence of dna of interest
		 * @param binWidth the width of the bins to be created
		 */
		public void createMatrix(DNAregion regionOfInterest,int binWidth) {
	
				//reads the start value of the sequence and end value of it
				long start=regionOfInterest.start;
				long end=regionOfInterest.end;
				//calculates the number of bins for the sequence
				int nrOfBins=(int) (end-start)/binWidth;
				//creates the bins and adding them to an array of type bin
				for(int i= (int)start;i<end; i=i+binWidth) {
				  DNAregion bin= new DNAregion(regionOfInterest.chromosome,i,i+binWidth);
				  bins.add(bin);
				}
		
				// creates the empty matrix with the values of 1 on diagonal
				int n=nrOfBins;
				CumulativeSimulatedMatrix = new double[n][n];
				ContactMatrix = new double[n][n];
				SimulatedMatrix= new boolean[n][n];
				NormalisedMatrix= new double[n][n];
				for(int i=0;i<nrOfBins;i++) { 
				  ContactMatrix[i][i]=1;
				  SimulatedMatrix[i][i]=true;
				  NormalisedMatrix[i][i]=0;
				  CumulativeSimulatedMatrix[i][i]=0;
				}
				
		 }
		
		
		/**
		 * Iterates through the simulated matrix and if the bins are interacting 
		 * the cumulative matrix will give the amount of time the bins were interacting
		 * @param newTime the time at which the simulation happened
		 */
		public void addValueToCumulativeSimulatedMatrix(double newTime) {
			double addTime= newTime-LastTimeSimulated;
			for(int i=0;i<bins.size();i++) {
				for(int j=0;j<bins.size();j++) {
					if(SimulatedMatrix[i][j]) {
						CumulativeSimulatedMatrix[i][j]+=addTime;
					}
				}
			}
		}
		
		/**
		 * reads the file of interest and adds the values from it to a local file
		 * the local file is read and the values for local arraylists for bins and score are added
		 * @param ContactMatrixLocation the path to the file containing the HiC interactions
		 * @throws Exception 
		*/
		public void addScoreToMatrix(DNAregion regionOfInterest, String ContactMatrixLocation) {
			 	//initialises every list
			 	ArrayList<String> localFile= new ArrayList<String>();
			 	ArrayList<DNAregion> binsX= new ArrayList<DNAregion>();
			 	ArrayList<DNAregion> binsY= new ArrayList<DNAregion>();		    
			 	ArrayList<Double> score = new ArrayList<Double>();
			 	
			 	//reads the interaction matrix file 
			 	String fileName= ContactMatrixLocation; 
			 	File file = new File(fileName);

			 	try {
				  Scanner inputStream= new Scanner(file);				  
				  while (inputStream.hasNextLine()) {
					  String data = inputStream.nextLine();
					  String[] values = data.split("\t");
					  for(String val:values) {
					     localFile.add(val);
					    
					  } 
			      }
				  inputStream.close();		 
			 	}  
			 	catch(FileNotFoundException e) {
				 e.printStackTrace();
			 	}
			 	//adds the values of the score from the local file to the arraylist score
			 	for(int i=6; i<localFile.size(); i=i+7) {
				     double value = Double.parseDouble(localFile.get(i));
					 score.add(value);
			 	}
			 	//local bins X and Y from file
			 	for(int i=1; i<localFile.size(); i=i+7) {
				  DNAregion bin= new DNAregion(localFile.get(0),Integer.parseInt(localFile.get(i)), Integer.parseInt(localFile.get(i+1)));
					binsX.add(bin);
			 	}
			 	for(int i=4; i<localFile.size(); i=i+7) {
				   DNAregion  bin = new DNAregion(localFile.get(3),Integer.parseInt(localFile.get(i)), Integer.parseInt(localFile.get(i+1)));
					binsY.add(bin);
			 	}
			 	boolean anyValue = false;
			
			 	//adds the values from the score for the specified region of interest into the matrix
			 	for(int g=0;g<score.size();g++) {
					  int indexBinX=getBinIndex(bins,binsX.get(g)); 
					  int indexBinY=getBinIndex(bins,binsY.get(g));
					  if(indexBinX>=0 && indexBinY>=0) { 
				    	ContactMatrix[indexBinX][indexBinY]=score.get(g);	
				    	anyValue = true;
				    	
					  }
			 	}
			 	//throws an error if there is a chromosome mismatch between the interaction matrix and the dna region
			 	if(!anyValue) {
					throw new IllegalArgumentException("Chromosome in the InteractionMatrix does not match the chromosome in the DNA sequence.");
			 	}	      
		 }
		 
         /**
          * the method looks into the array and compares each object in the array with the object to be found
          * returns the index in the array at which the bin has the same value
          * @param bins array of type Bin  
          * @param subject object of class Bin
          * @return the bin position in the list based on its details
          */
		 private int getBinIndex(ArrayList<DNAregion> bins,DNAregion subject) {
			 	for(int i=0; i<bins.size(); i++) {
			      if(bins.get(i).equals(subject)) {
				     return i;
				  }
			 	}
			 	return -1;
		 }
		 
		 /**
		  * Method that searches the bin in which a certain position is in
		  * @param position the position of the object of interest
		  * @param chromosome the input param has to match the chromosome of the bins from the list
		  * @return the index of the current bin
		  */
		 public int getCurrentBinIndex(int position, String chromosome) {
				int local=-1;
				for(int i = 0; i<bins.size(); i++) {
					if(chromosome.equals(bins.get(i).chromosome) && position>=bins.get(i).start && position<=bins.get(i).end) {
					 local= i;

					}
				}
				return local;		
		 }
		 
		 /**
		  * gets the bin's details based on its index
		  * @param index the index of the bin
		  * @return the details of the bin
		  */
		 public DNAregion getBin(int index) {
			 	return bins.get(index);
		 }
		
		 /**
		  * takes the input parameter and searches against the existing arraylist of bins to get the interacting bins
		  * @param bin the no of the bin of interest
		  * @return a list of indexes of bins
		  */
		 public ArrayList<Integer> getInteractingBins(int bin){
			 
			 	ArrayList<Integer> local = new ArrayList<Integer>();
			 	for(int i =0; i<bins.size(); i++) {
					 if(SimulatedMatrix[bin][i]) {
						local.add(i);
					 }
			 	}
			 	return local;
		 }
		 
		 /**
		  * gives the size of the bins list
		  * @return integer = the size of bins list
		  */
		 public int getBinSize() {
			 	return bins.size();
		 }
		 		 
		 /**
		  * prints out details regarding 2 bins and the score between them as well as the bins ID
		  */
		 public void find(int BinX,int BinY) {
			 	System.out.println("BinX: "+BinX+ "BinY: "+BinY +"Score: " +ContactMatrix[BinX][BinY] +"Bins ID: " +bins.get(BinX) +" "+ bins.get(BinY) );
		 }
		 
		 /**
		  * Accessor method for the arraylist 'bins';
		  * @return a list of type ArrayList<DNAregion>
		  */
		 public ArrayList<DNAregion> getBinsList() {
			 	return bins;
		 }
		 
		 /**
		  * Method that takes the highest score from the score list
		  * and divides every value from the interaction matrix by the highest score to get the normalised values between 0 and 1
		  */
		 private void normaliseMatrix() {
			 	double highestScore = 0;
			 	//get the highest score
			 	for(int i=0;i<bins.size();i++) {
				 for(int j=0;j<bins.size();j++) {
					 if(ContactMatrix[i][j]>highestScore) {
						 highestScore=ContactMatrix[i][j];
					 }
				 }	 
			 	}
			 	//generate normalised matrix
			 	for(int i=0;i<bins.size();i++) {
				 for(int j=0;j<bins.size();j++) {
					 NormalisedMatrix[i][j]=ContactMatrix[i][j]/highestScore;
				 }
			 	}
		 }
		 
		 /**
		  * Simulates the Interaction Matrix
		  * Can be accessed and simulated from other classes to have a dynamic simulation matrix
		  * @param generator
		  */
		 public void simulateMatrix(Random generator, double time) {		
			 LastTimeSimulated=time;
				for(int i=0;i<bins.size();i++) {
					for(int j=0; j<bins.size();j++) {
						double probability=Utils.generateNextDouble(generator, 0, 1);
						SimulatedMatrix[i][j]=(NormalisedMatrix[i][j]>probability);
					}
				}
		 }
		 
		 /**
		  * method that returns the cumulative simulation time between two bins
		  * @param t bin X
		  * @param b bin Y
		  * @return 
		  */
		 public double CumulativeResult(int t, int b) {
			 return CumulativeSimulatedMatrix[t][b];
		 }
		 
		 /**
		  * method that generates a random number between the values of the start and end of a bin
		  */
		 public int radomBinPosition(Random generator,int randomBinID) {
		  int newPosition = (generator.nextInt((int)(bins.get(randomBinID).end-bins.get(randomBinID).start)+1)+(int)bins.get(randomBinID).start)-(int)bins.get(0).start;
		 return newPosition;
		 }
		 
		 /**
		  * 
		  * @return
		  */
		 public String headerToString() {
			 String str = "\"BIN_X\", \"BIN_Y\", \"HIC_SCORE\", \"CUMULATIVE_SIMULATION\"";
				return str;
		 }
		 
		 /**
		  * a string with the details about the bins interacting, their interaction score and the cumulative time they are interacting
		  * @param x bin X
		  * @param y bin Y
		  * @return
		  */
		 public String toString(int x, int y){
				//"\"BIN_X\", \"BIN_Y\",
				String str="\""+ bins.get(x);
				//\"BIN_Y\", \"HIC_SCORE\", \"CUMULATIVE_SIMULATION\"
				str+="\","+bins.get(y)+", "+ContactMatrix[x][y]+", "+CumulativeSimulatedMatrix[x][y];
				return str;
			}
}
