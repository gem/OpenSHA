package org.opensha.sha.imr.attenRelImpl;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;

import org.opensha.commons.data.NamedObjectAPI;
import org.opensha.commons.data.Site;
import org.opensha.commons.param.DoubleConstraint;
import org.opensha.commons.param.DoubleDiscreteConstraint;
import org.opensha.commons.param.StringConstraint;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeListener;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.attenRelImpl.constants.ToroEtAl2002Constants;
import org.opensha.sha.imr.param.EqkRuptureParams.FaultTypeParam;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.EqkRuptureParams.RakeParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.DampingParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGV_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.ComponentParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceJBParameter;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

/**
 * <b>Title:</b> Toro_2002_AttenRel
 * <p>
 * 
 * <b>Description:</b> This implements the updated GMPE developed by Toro et al
 * (1997) with the distance term proposed by Toro (2002)
 * (www.riskeng.com/PDF/atten_toro_extended.pdf)
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
 * <LI>componentParam - average horizontal , average horizontal (GMRoti50)
 * <LI>stdDevTypeParam - total, none
 * </UL>
 * <p>
 * 
 * @author l.danciu
 * @created July, 2011
 * @version 1.0
 */

public class ToroEtAl_2002_AttenRel extends AttenuationRelationship implements
		ScalarIntensityMeasureRelationshipAPI, NamedObjectAPI,
		ParameterChangeListener {

	// Debugging stuff
	private static String C = "ToroEtAl_2002_AttenRel";
	private static boolean D = false;
	public static String SHORT_NAME = "ToroEtAl2002";
	private static long serialVersionUID = 1234567890987654353L;

	// Name of IMR
	public static String NAME = "Toro et al. (2002)";

	/** Period index. */
	protected int iper;

	/** Moment magnitude. */
	protected double mag;

	/** Rupture distance. */
	protected double rJB;

	/** Standard deviation type. */
	protected String stdDevType;

	/** Map period-value/period-index. */
	protected HashMap<Double, Integer> indexFromPerHashMap;

	/** For issuing warnings. */
	private transient ParameterChangeWarningListener warningListener = null;

	/**
	 * Construct attenuation relationship. Initialize parameters and parameter
	 * lists.
	 */
	public ToroEtAl_2002_AttenRel(ParameterChangeWarningListener warningListener) {

		super();

		this.warningListener = warningListener;

		initSupportedIntensityMeasureParams();
		indexFromPerHashMap = new HashMap<Double, Integer>();
		for (int i = 1; i < ToroEtAl2002Constants.PERIOD.length; i++) {
			indexFromPerHashMap.put(
					new Double(ToroEtAl2002Constants.PERIOD[i]),
					new Integer(i));
		}

		initEqkRuptureParams();
		initPropagationEffectParams();
		initSiteParams();
		initOtherParams();
		initIndependentParamLists();
		initParameterEventListeners();
	}

	/**
	 * Creates the two supported IM parameters (PGA and SA), as well as the
	 * independenParameters of SA (periodParam and dampingParam) and adds them
	 * to the supportedIMParams list. Makes the parameters non-editable.
	 */
	protected void initSupportedIntensityMeasureParams() {

		// set supported periods for spectral acceleration
		DoubleDiscreteConstraint periodConstraint = new DoubleDiscreteConstraint();
		for (int i = 1; i < ToroEtAl2002Constants.PERIOD.length; i++) {
			periodConstraint.addDouble(new Double(
					ToroEtAl2002Constants.PERIOD[i]));
		}
		periodConstraint.setNonEditable();
		// set period param (default is 1s, which is provided by Toro et al 2002
		// GMPE)
		saPeriodParam = new PeriodParam(periodConstraint);

		// set damping parameter. Empty constructor set damping
		// factor to 5 % (which is the one provided by Toro et 2002 GMPE)
		saDampingParam = new DampingParam();

		// initialize spectral acceleration parameter (units: g)
		saParam = new SA_Param(saPeriodParam, saDampingParam);
		saParam.setNonEditable();

		// initialize peak ground acceleration parameter (units: g)
		pgaParam = new PGA_Param();
		pgaParam.setNonEditable();

		// add the warning listeners
		saParam.addParameterChangeWarningListener(warningListener);
		pgaParam.addParameterChangeWarningListener(warningListener);

		// put parameters in the supportedIMParams list
		supportedIMParams.clear();
		supportedIMParams.addParameter(saParam);
		supportedIMParams.addParameter(pgaParam);

	}

	/**
	 * Initialize earthquake rupture parameter (moment magnitude, rake) and add
	 * to eqkRuptureParams list. Makes the parameters non-editable.
	 */
	protected void initEqkRuptureParams() {

		// moment magnitude (default 5.0)
		magParam = new MagParam(ToroEtAl2002Constants.MAG_WARN_MIN,
				ToroEtAl2002Constants.MAG_WARN_MAX);
		eqkRuptureParams.clear();
		eqkRuptureParams.addParameter(magParam);

	}

	/**
	 * Initialize site parameters. No site parameters are defined so the site
	 * parameter list is made empty.
	 */
	protected void initSiteParams() {
		siteParams.clear();
	}

	/**
	 * Initialize Propagation Effect parameters (JB distance) and adds them to
	 * the propagationEffectParams list. Makes the parameters non-editable.
	 */
	protected void initPropagationEffectParams() {

		distanceJBParam = new DistanceJBParameter(
				ToroEtAl2002Constants.DISTANCE_JB_WARN_MIN);
		distanceJBParam.addParameterChangeWarningListener(warningListener);
		DoubleConstraint warn = new DoubleConstraint(new Double(0.00),
				ToroEtAl2002Constants.DISTANCE_JB_WARN_MAX);
		warn.setNonEditable();
		distanceJBParam.setWarningConstraint(warn);
		distanceJBParam.setNonEditable();

		propagationEffectParams.addParameter(distanceJBParam);
	}

	/**
	 * Initialize other Parameters (standard deviation type, component, sigma
	 * truncation type, sigma truncation level).
	 */
	protected void initOtherParams() {

		sigmaTruncTypeParam = new SigmaTruncTypeParam();
		sigmaTruncLevelParam = new SigmaTruncLevelParam();

		// stdDevType Parameter
		StringConstraint stdDevTypeConstraint = new StringConstraint();
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_TOTAL);
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_NONE);
		stdDevTypeConstraint.setNonEditable();
		stdDevTypeParam = new StdDevTypeParam(stdDevTypeConstraint);

		// the component Parameter
		StringConstraint constraint = new StringConstraint();
		constraint.addString(ComponentParam.COMPONENT_AVE_HORZ);
		constraint.addString(ComponentParam.COMPONENT_GMRotI50);
		constraint.setNonEditable();
		componentParam = new ComponentParam(constraint,
				ComponentParam.COMPONENT_AVE_HORZ);

		// add these to the list
		otherParams.clear();
		otherParams.addParameter(sigmaTruncTypeParam);
		otherParams.addParameter(sigmaTruncLevelParam);
		otherParams.addParameter(stdDevTypeParam);
		otherParams.addParameter(componentParam);
	}

	/**
	 * This creates the lists of independent parameters that the various
	 * dependent parameters (mean, standard deviation, exceedance probability,
	 * and IML at exceedance probability) depend upon. NOTE: these lists do not
	 * include anything about the intensity-measure parameters or any of their
	 * internal independentParamaters.
	 */
	protected void initIndependentParamLists() {

		// params that the mean depends upon
		meanIndependentParams.clear();
		meanIndependentParams.addParameter(magParam);
		meanIndependentParams.addParameter(distanceJBParam);

		// params that the stdDev depends upon
		stdDevIndependentParams.clear();
		stdDevIndependentParams.addParameter(stdDevTypeParam);

		// params that the exceed. prob. depends upon
		exceedProbIndependentParams.clear();
		exceedProbIndependentParams.addParameterList(meanIndependentParams);
		exceedProbIndependentParams.addParameter(stdDevTypeParam);
		exceedProbIndependentParams.addParameter(sigmaTruncTypeParam);
		exceedProbIndependentParams.addParameter(sigmaTruncLevelParam);

		// params that the IML at exceed. prob. depends upon
		imlAtExceedProbIndependentParams
				.addParameterList(exceedProbIndependentParams);
		imlAtExceedProbIndependentParams.addParameter(exceedProbParam);
	}

	/**
	 * This sets the eqkRupture related parameters (moment magnitude, tectonic
	 * region type, focal depth) based on the eqkRupture passed in. The
	 * internally held eqkRupture object is also set as that passed in. Warning
	 * constrains on magnitude and focal depth are ignored.
	 */
	public void setEqkRupture(EqkRupture eqkRupture) {

		magParam.setValueIgnoreWarning(new Double(eqkRupture.getMag()));
		this.eqkRupture = eqkRupture;
		setPropagationEffectParams();
	}

	/**
	 * Sets the internally held Site object as that passed in.
	 */
	public void setSite(Site site) {
		this.site = site;
		setPropagationEffectParams();
	}

	/**
	 * This calculates the Rupture Distance propagation effect parameter based
	 * on the current site and eqkRupture.
	 * <P>
	 */
	protected void setPropagationEffectParams() {

		if ((this.site != null) && (this.eqkRupture != null)) {
			distanceJBParam.setValue(eqkRupture, site);
		}
	}

	/**
	 * Set period index.
	 */
	protected void setPeriodIndex() {
		if (im.getName().equalsIgnoreCase(PGA_Param.NAME)) {
			iper = 0;
		} else {
			iper = ((Integer) indexFromPerHashMap.get(saPeriodParam.getValue()))
					.intValue();
		}
	}

	/**
	 * Compute mean.
	 * 
	 */
	public double getMean() {
		if (rJB > USER_MAX_DISTANCE) {
			return VERY_SMALL_MEAN;
		} else {
			setPeriodIndex();
			return getMean(iper, mag, rJB);
		}
	}

	public double getStdDev() {
		setPeriodIndex();
		return getStdDev(iper, mag, rJB, stdDevType);
	}

	/**
	 * This computes the mean ln(Y)
	 * 
	 * @param iper
	 * @param rJB
	 * @param mag
	 * @param vs30
	 * @param rake
	 */
	public double getMean(int iper, double mag, double rJB) {
		/**
		 * This is to avoid very small values for rJB
		 * 
		 * */
		if (rJB < 1e-3) {
			rJB = 1;
		}

		double magDiff = mag - 6.0;

		double rM = Math.sqrt(rJB * rJB + ToroEtAl2002Constants.c7[iper]
				* ToroEtAl2002Constants.c7[iper]
				* Math.pow(Math.exp(-1.25 + 0.227 * mag), 2));

		double f1 = ToroEtAl2002Constants.c1[iper]
				+ ToroEtAl2002Constants.c2[iper] * magDiff
				+ ToroEtAl2002Constants.c3[iper] * magDiff * magDiff;

		double f2 = ToroEtAl2002Constants.c4[iper] * Math.log(rM);

		double f3 = (ToroEtAl2002Constants.c5[iper] - ToroEtAl2002Constants.c4[iper])
				* Math.max(Math.log(rM / 100), 0);

		double f4 = ToroEtAl2002Constants.c6[iper] * rM;

		double lnY = f1 - f2 - f3 - f4;
		
		// correction for 3 s
		if  (ToroEtAl2002Constants.PERIOD[iper] == 3.00) {
			lnY /= ToroEtAl2002Constants.T2sec_TO_T3sec_DECAYFACTOR;
		}
		// correction for 4 s
		if  (ToroEtAl2002Constants.PERIOD[iper] == 4.00) {
			lnY /= ToroEtAl2002Constants.T2sec_TO_T4sec_DECAYFACTOR;
		}

		return lnY;

	}

	public double getStdDev(int iper, double mag, double rJB, String stdDevType) {
		double sigmaaM = Double.NaN;
		double sigmaaR = Double.NaN;
		double sigmae = Double.NaN;
		double sigmatot = Double.NaN;

		if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_NONE))
			return 0;
		else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_TOTAL)) {
			if (mag < 5.0) {
				sigmaaM = ToroEtAl2002Constants.m50[iper];
			} else if ((mag >= 5.0) && (mag <= 5.5)) {
				sigmaaM = ToroEtAl2002Constants.m50[iper]
						+ (ToroEtAl2002Constants.m55[iper] - ToroEtAl2002Constants.m50[iper])
						/ (5.5 - 5.0) * (mag - 5.0);
			} else if ((mag > 5.5) && (mag < 8.0)) {
				sigmaaM = ToroEtAl2002Constants.m55[iper]
						+ (ToroEtAl2002Constants.m80[iper] - ToroEtAl2002Constants.m55[iper])
						/ (8.0 - 5.5) * (mag - 5.5);
			} else {
				sigmaaM = ToroEtAl2002Constants.m80[iper];
			}

			if (rJB < 5.0) {
				sigmaaR = ToroEtAl2002Constants.r05[iper];
			} else if ((rJB >= 5.0) && (rJB <= 20.0)) {
				sigmaaR = ToroEtAl2002Constants.r05[iper]
						+ (ToroEtAl2002Constants.r20[iper] - ToroEtAl2002Constants.r05[iper])
						/ (20.0 - 5.0) * (rJB - 5.0);
			} else {
				sigmaaR = ToroEtAl2002Constants.r20[iper];
			}

			if (ToroEtAl2002Constants.PERIOD[iper] >= 2.0) {
				sigmae = 0.34 + 0.06 * (mag - 6.0);
			} else {
				sigmae = 0.36 + 0.07 * (mag - 6.0);
			}
			double sigmaatot = Math.sqrt(sigmaaM * sigmaaM + sigmaaR * sigmaaR);

			sigmatot = (Math.sqrt(sigmaatot * sigmaatot + sigmae * sigmae));

			return sigmatot;
		} else
			throw new RuntimeException("Standard deviation type: " + stdDevType
					+ " not recognized");
	}

	/**
	 * Allows the user to set the default parameter values for the selected
	 * Attenuation Relationship.
	 */
	public void setParamDefaults() {

		magParam.setValueAsDefault();
		distanceJBParam.setValueAsDefault();
		saPeriodParam.setValueAsDefault();
		saDampingParam.setValueAsDefault();
		saParam.setValueAsDefault();
		pgaParam.setValueAsDefault();
		stdDevTypeParam.setValueAsDefault();
		sigmaTruncTypeParam.setValueAsDefault();
		sigmaTruncLevelParam.setValueAsDefault();
		componentParam.setValueAsDefault();
	}

	/**
	 * This listens for parameter changes and updates the primitive parameters
	 * accordingly
	 */
	public void parameterChange(ParameterChangeEvent e) {

		String pName = e.getParameterName();
		Object val = e.getNewValue();

		if (pName.equals(MagParam.NAME)) {
			mag = ((Double) val).doubleValue();
		} else if (pName.equals(DistanceJBParameter.NAME)) {
			rJB = ((Double) val).doubleValue();
		} else if (pName.equals(StdDevTypeParam.NAME)) {
			stdDevType = (String) val;
		}
	}

	/**
	 * Allows to reset the change listeners on the parameters
	 */
	public void resetParameterEventListeners() {
		magParam.removeParameterChangeListener(this);
		distanceJBParam.removeParameterChangeListener(this);
		stdDevTypeParam.removeParameterChangeListener(this);
		saPeriodParam.removeParameterChangeListener(this);
		this.initParameterEventListeners();
	}

	/**
	 * Adds the parameter change listeners. This allows to listen to when-ever
	 * the parameter is changed.
	 */
	protected void initParameterEventListeners() {

		magParam.addParameterChangeListener(this);
		distanceJBParam.addParameterChangeListener(this);
		stdDevTypeParam.addParameterChangeListener(this);
		saPeriodParam.addParameterChangeListener(this);
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

	/**
	 * Returns the URL of the AttenuationRelationship documentation. Currently
	 * returns null because no URL has been created yet.
	 */
	public URL getAttenuationRelationshipURL() throws MalformedURLException {
		return null;
	}
}
