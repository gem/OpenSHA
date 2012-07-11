package org.opensha.sha.imr.attenRelImpl;

import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.attenRelImpl.constants.AdjustFactorsSHARE;
import org.opensha.sha.imr.param.EqkRuptureParams.RakeParam;

/**
 * <b>Title:</b> Campbell2003share_AttenRel
 * <p>
 * 
 * <b>Description:</b> Class implementing GMPE described in: Prediction of
 * Strong Ground Motion Using the Hybrid Empirical Method and Its Use in the
 * Development of Ground-Motion (Attenuation) Relations in Eastern North America
 * (June 2003, BSSA, vol 93, no 3, pp 1012-1033). Modifications of the equations
 * according to the ERRATUM (July 2004) are also included.
 * <p>
 * The GMPE is adjusted to account for style of faulting and a default rock soil
 * (Vs30 >=800m/sec). The adjustment coefficients were proposed by S. Drouet
 * [2010] - internal SHARE WP4 report;
 * <p>
 * <UL>
 * Supported Intensity-Measure Parameters:
 * <LI>PGA - Peak Ground Acceleration
 * <LI>saParam - Response Spectral Acceleration
 * <LI>
 * </UL>
 * <p>
 * Other Independent Parameters:
 * <p>
 * <UL>
 * <LI>magParam - moment magnitude
 * <LI>distanceRupParam - closest distance to rupture
 * <LI>rakeParam - rake angle
 * <LI>componentParam - average horizontal, average horizontal (GMRoti50)
 * <LI>stdDevTypeParam - none, total
 * </UL>
 * <p>
 * 
 * @author l.danciu
 * @created October, 2010 - updated July 2011
 * @version 1.01
 */
public class Campbell_2003_SHARE_AttenRel extends Campbell_2003_AttenRel {

	public final static String NAME = "Campbell 2003 (SHARE)";

	public final static String SHORT_NAME = "Campbell2003(SHARE)";

	private static final long serialVersionUID = 1L;

	private double rake;

	public void setIper(int iper) {
	    this.iper = iper;
	}

	public Campbell_2003_SHARE_AttenRel(
			ParameterChangeWarningListener warningListener) {
		super(warningListener);
	}

	/**
	 * Call superclass method and initialize rake parameter from earthquake
	 * rupture.
	 */
	protected void initEqkRuptureParams() {
		super.initEqkRuptureParams();
		rakeParam = new RakeParam();
		eqkRuptureParams.addParameter(rakeParam);
	}

	/**
	 * Call superclass method and add rake parameter to the list of parameters
	 * the mean calculation depends upon.
	 */
	protected void initIndependentParamLists() {
		super.initIndependentParamLists();
		meanIndependentParams.addParameter(rakeParam);
	}

	/**
	 * This sets the eqkRupture related parameters (magnitude and rake angle)
	 * based on the eqkRupture passed in. The internally held eqkRupture object
	 * is also set as that passed in.
	 */
	public void setEqkRupture(EqkRupture eqkRupture) {
		magParam.setValueIgnoreWarning(new Double(eqkRupture.getMag()));
		if (!Double.isNaN(eqkRupture.getAveRake())) {
			rakeParam.setValue(eqkRupture.getAveRake());
		}
		this.eqkRupture = eqkRupture;
		setPropagationEffectParams();
	}

	/**
	 * Compute ground motion distribution mean value (in ln space).
	 */
	public double getMean() {
		if (rRup > USER_MAX_DISTANCE) {
			return VERY_SMALL_MEAN;
		} else {
			setPeriodIndex();
			return getMean(iper, mag, rRup, rake);
		}
	}

	/**
	 * Compute adjusted ground motion distribution mean value (in ln space)
	 */
	public double getMean(int iper, double mag, double rRup, double rake) {
		double meanOriginal = super.getMean(iper, mag, rRup);
		double[] f = computeStyleOfFaultingTerm(iper, rake);
		return Math.log(Math.exp(meanOriginal) * f[2]
				* AdjustFactorsSHARE.AFrock_CAMPBELL2003[iper]);
	}

	/**
	 * Compute style-of-faulting adjustment
	 **/
	public double[] computeStyleOfFaultingTerm(final int iper, final double rake) {
		double[] f = new double[3];
		if (rake > AdjustFactorsSHARE.FLT_TYPE_NORMAL_RAKE_LOWER
				&& rake <= AdjustFactorsSHARE.FLT_TYPE_NORMAL_RAKE_UPPER) {
			f[0] = 1.0;
			f[1] = 0.0;
			f[2] = f[0]
					* Math.pow(AdjustFactorsSHARE.Frss_CAMPBELL2003[iper],
							(1 - AdjustFactorsSHARE.pR))
					* Math.pow(AdjustFactorsSHARE.Fnss, -AdjustFactorsSHARE.pN);
		} else if (rake > AdjustFactorsSHARE.FLT_TYPE_REVERSE_RAKE_LOWER
				&& rake <= AdjustFactorsSHARE.FLT_TYPE_REVERSE_RAKE_UPPER) {
			f[0] = 0.0;
			f[1] = 1.0;
			f[2] = f[1]
					* Math.pow(AdjustFactorsSHARE.Frss_CAMPBELL2003[iper],
							-AdjustFactorsSHARE.pR)
					* Math.pow(AdjustFactorsSHARE.Fnss,
							(1 - AdjustFactorsSHARE.pN));
		} else {
			f[0] = 0.0;
			f[1] = 0.0;
			f[2] = Math.pow(AdjustFactorsSHARE.Frss_CAMPBELL2003[iper],
					-AdjustFactorsSHARE.pR)
					* Math.pow(AdjustFactorsSHARE.Fnss, -AdjustFactorsSHARE.pN);
		}
		return f;
	}

	public double getStdDev(){
		setPeriodIndex();
	    return getStdDev(iper, mag, stdDevType);
	}

	public double getStdDev(int iper, double mag, String stdDevType){
		double std = super.getStdDev(iper, mag, stdDevType);
		return std * AdjustFactorsSHARE.sig_AFrock_CAMBPELL2003[iper];
	}

	/**
	 * Allows the user to set the default parameter values for the selected
	 * Attenuation Relationship.
	 */
	public void setParamDefaults() {
		super.setParamDefaults();
		rakeParam.setValueAsDefault();
	}

	/**
	 * This listens for parameter changes and updates the primitive parameters
	 * accordingly
	 */
	public void parameterChange(ParameterChangeEvent e) {

		super.parameterChange(e);

		String pName = e.getParameterName();
		Object val = e.getNewValue();

		if (pName.equals(RakeParam.NAME)) {
			rake = ((Double) val).doubleValue();
		}
	}

	/**
	 * Allows to reset the change listeners on the parameters
	 */
	public void resetParameterEventListeners() {
		super.resetParameterEventListeners();
		rakeParam.removeParameterChangeListener(this);
		rakeParam.addParameterChangeListener(this);
	}

	/**
	 * Get the name of this IMR.
	 */
	public String getName() {
		return NAME;
	}

	/**
	 * Returns the Short Name of each AttenuationRelationship
	 * 
	 */
	public String getShortName() {
		return SHORT_NAME;
	}

}
