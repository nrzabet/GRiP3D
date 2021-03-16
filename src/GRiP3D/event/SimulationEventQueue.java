package event;

import java.io.Serializable;

import environment.Cell;

/**
 * class that holds all the simulation events
 * @author adumita@essex.ac.uk
 *
 */
public class SimulationEventQueue implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	protected double simulationPropensity;// the propensity that the matrix will be simulated again;
	
	private SimulationEvent simulationEvent;
	
	public SimulationEventQueue(Cell n){
				this.simulationPropensity=n.ip.PROPORTION_TIME.value;
				this.simulationEvent = null;
				//simulationEvent.time= generateExponentialDistribution(simulationPropensity,n);
	}
	 
	public double generateExponentialDistribution(double mean,Cell n) {
		return Math.log(1-n.randomGenerator.nextDouble())/(-mean);
	}
	
	/**
	 * returns the next Simulation event
	 * @return
	 */
	public SimulationEvent peek(){
		return simulationEvent;
	}
	
	/**
	 * returns the next Simulation event and removes it from the list
	 * @return
	 */
	public SimulationEvent pop(){
		SimulationEvent pe=simulationEvent;
		simulationEvent = null;
		return pe;
	}
	
	/**
	 * replaces the current Simulation event with a new one
	 * @param re the new event
	 */
	public void add(SimulationEvent pe){
		simulationEvent=pe;
		//System.out.println(re);
	}
	
	
	/**
	 * returns true if the Simulation event is null
	 * @return
	 */
	public boolean isEmpty(){
		return simulationEvent==null;
	}
	
	
	
	/**
	 * deletes current protein binding event
	 */
	public void clear(){
		this.simulationEvent = null;
	}

}
