package event;

/**
 *class that describes a simulation event. it is an instantiation of Event
 * @author adumita@essex.ac.uk
 *
 */
public class SimulationEvent extends Event {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public boolean isSimulationEvent;
	
	public SimulationEvent(double time, int nextAction, boolean isSimulationEvent) {
		super(time,nextAction);
		this.isSimulationEvent = isSimulationEvent;
	}
	
	
	/**
	 * generates the description string of current event
	 */
	public String toString(){
		String stateStr=""+time+": ";
		stateStr+=" through an event of type "+nextAction;
		return stateStr;
	}
	/**
	 * compares whether this event equals the one supplied as an argument
	 * @param pe
	 * @return
	 */
	public boolean isEqualTo(SimulationEvent pe){
		return this.time == pe.time;
	}
}
