package org.opensha.sha.imr.attenRelImpl;

import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.attenRelImpl.constants.AdjustFactorsSHARE;
import org.opensha.sha.imr.attenRelImpl.constants.AdjustFactorsSHARE;
import org.opensha.sha.imr.attenRelImpl.constants.ToroEtAl2002Constants;
import org.opensha.sha.imr.param.EqkRuptureParams.RakeParam;

/**
 * <b>Title:</b> ToroEtAl_2002_SHARE_AttenRel
 * <p>
 * 
 * <b>Description:</b> This class implements the updated GMPE developed by Toro
 * et al - 2002, (www.riskeng.com/PDF/atten_toro_extended.pdf): MODIFICATION OF
 * THE TORO ET AL. (1997) ATTENUATION EQUATIONS FOR LARGE MAGNITUDES AND SHORT
 * DISTANCES.
 * 
 * The GMPE is adjusted to account the style of faulting and a default rock soil
 * (Vs30 >=800m/sec) The adjustment coefficients were proposed by S. Drouet
 * [2010]; Supported period values (s). Period 0.5s was obtained as a linear
 * interpolation between 0.4 and 1.0s; Warning: The coefficients for periods
 * 3.00 and 4.00sec were obtained as a function of SA(2sec) and ratios between
 * SA(3)/SA(4sec) and SA(2)/SA(4sec) Disclaimer: The adjustment of the SA(3sec)
 * and SA(4sec) are obtained in the framework of SHARE project.
 * 
 * Supported Intensity-Measure Parameters:
 * <p>
 * <UL>
 * <LI>saParam - Response Spectral Acceleration
 * <LI>PGA - Peak Ground Acceleration
 * </UL>
 * <p>
 * Other Independent Parameters:
 * <p>
 * <UL>
 * <LI>magParam - moment magnitude
 * <LI>distanceJBParam - JB distance
 * <LI>componentParam - average horizontal, average horizontal (GMRoti50)
 * <LI>rakeParam - rake angle
 * <LI>stdDevTypeParam - total, none
 * </UL>
 * <p>
 * 
 * @author l.danciu
 * @created July, 2011
 * @version 1.0
 */
public class ToroEtAl_2002_SHARE_AttenRel extends ToroEtAl_2002_AttenRel {

	private static final long serialVersionUID = 1L;

	public final static String NAME = "Toro et al. (2002) (SHARE)";

	public final static String SHORT_NAME = "ToroEtAl2002(SHARE)";

	private double rake;

	/**
	 * Calls the {@link ToroEtAl_2002_AttenRel} constructor.
	 */
	public ToroEtAl_2002_SHARE_AttenRel(
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
		if (rJB > USER_MAX_DISTANCE) {
			return VERY_SMALL_MEAN;
		} else {
			setPeriodIndex();
			return getMean(iper, mag, rJB, rake);
		}
	}

	/**
	 * Compute adjusted mean value taking into account the faulting style and
	 * new reference value for rock (vs30>=800).
	 */
	public double getMean(int iper, double mag, double rJB, double rake) {
		double meanOriginal = super.getMean(iper, mag, rJB);
		double meanAdj = Double.NaN;
		double[] f = computeStyleOfFaultingTerm(iper, rake);
		meanAdj = Math.exp(meanOriginal) * f[2]
				* AdjustFactorsSHARE.AFrock_TORO2002[iper];
		return Math.log(meanAdj);
	}

	/**
	 * Compute style-of-faulting adjustment
	 **/
	public double[] computeStyleOfFaultingTerm(final int iper, final double rake) {
		double[] f = new double[3];
		if (rake > ToroEtAl2002Constants.FLT_TYPE_NORMAL_RAKE_LOWER
				&& rake <= ToroEtAl2002Constants.FLT_TYPE_NORMAL_RAKE_UPPER) {
			f[0] = 1.0;
			f[1] = 0.0;
			f[2] = f[0]
					* Math.pow(AdjustFactorsSHARE.Frss[iper],
							(1 - AdjustFactorsSHARE.pR))
					* Math.pow(AdjustFactorsSHARE.Fnss,
							-AdjustFactorsSHARE.pN);
		} else if (rake > ToroEtAl2002Constants.FLT_TYPE_REVERSE_RAKE_LOWER
				&& rake <= ToroEtAl2002Constants.FLT_TYPE_REVERSE_RAKE_UPPER) {
			f[0] = 0.0;
			f[1] = 1.0;
			f[2] = f[1]
					* Math.pow(AdjustFactorsSHARE.Frss[iper],
							-AdjustFactorsSHARE.pR)
					* Math.pow(AdjustFactorsSHARE.Fnss,
							(1 - AdjustFactorsSHARE.pN));
		} else {
			f[0] = 0.0;
			f[1] = 0.0;
			f[2] = Math.pow(AdjustFactorsSHARE.Frss[iper],
					-AdjustFactorsSHARE.pR)
					* Math.pow(AdjustFactorsSHARE.Fnss,
							-AdjustFactorsSHARE.pN);
		}
		return f;
	}

	public double getStdDev() {
		return getStdDev(iper, mag, rJB, stdDevType);
	}

	public double getStdDev(int iper, double mag, double rJB, String stdDevType) {
		return super.getStdDev(iper, mag, rJB, stdDevType)
				* AdjustFactorsSHARE.sig_AFrock_TORO2002[iper];
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
