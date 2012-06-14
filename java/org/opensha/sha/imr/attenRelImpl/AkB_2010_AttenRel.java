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
import org.opensha.sha.imr.attenRelImpl.constants.AkB2010Constants;
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
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

/**
 * <b>Title:</b> AkB_2010_AttenRel
 * <p>
 * 
 * <b>Description:</b> Class implementing attenuation relationship described in:
 * "Empirical Equations for the Prediction of PGA , PGV, and Spectral
 * Accelerations in Europe, the Mediterranean Region, and the Middle East",
 * Sinnan Akkar, Julian J. Bommer, Seismological Research Letters Vol. 81, No 2,
 * pp 195-206, March-April 2010.
 * <p>
 * 
 * Supported Intensity-Measure Parameters:
 * <p>
 * <UL>
 * <LI>pgaParam - Peak Ground Acceleration (g)
 * <LI>pgvParam - Peak Ground Velocity (cm/s)
 * <LI>saParam - Response Spectral Acceleration (g)
 * </UL>
 * <p>
 * Other Independent Parameters:
 * <p>
 * <UL>
 * <LI>magParam - moment magnitude
 * <LI>rakeParam - rake angle. Used to establish if event is normal (-150 < rake
 * < -30) or reverse ((30 < rake < 150).
 * <LI>distanceJBParam - Joyner-Boore distance
 * <LI>vs30Param - shear wave velocity (m/s) averaged over the top 30 m of the
 * soil profile; The model assumes the following classification: vs30 < 360 m/s
 * -> soft soil, 360 <= vs30 <= 750 -> stiff soil, vs30 > 750 -> rock.
 * <LI>componentParam - average horizontal (geometrical mean), GMRoti50
 * (original model assumes only average horizontal, GMRorti50 added assuming
 * equivalence to average horizontal according to SHARE Report D4.2: Adjustment
 * of GMPEs, S. Drouet, F. Cotton, C. Beauval)
 * <LI>stdDevTypeParam - total, inter-event, intra-event, none
 * </UL>
 * <p>
 * 
 * <p>
 * Verification - Class tested (see {@link AkB_2010_test}) using tables provided
 * by the original authors
 * </p>
 * 
 ** 
 * @author L. Danciu, D. Monelli
 * @created August 20, 2010
 * @version 1.0
 */

public class AkB_2010_AttenRel extends AttenuationRelationship implements
		ScalarIntensityMeasureRelationshipAPI, NamedObjectAPI,
		ParameterChangeListener {

	/** Short name. */
	public final static String SHORT_NAME = "AkB2010";

	/** Full name. */
	public final static String NAME = "Akkar & Bommer 2010";

	/** Version number. */
	private static final long serialVersionUID = 1234567890987654353L;

	/** Period index. */
	private int iper;

	/** Moment magnitude. */
	private double mag;

	/** rake angle. */
	private double rake;

	/** Joyner and Boore distance. */
	private double rJB;

	/** Vs 30. */
	private double vs30;

	/** Standard deviation type. */
	private String stdDevType;

	/** Map period-value/period-index. */
	private HashMap<Double, Integer> indexFromPerHashMap;

	/** For issuing warnings. */
	private transient ParameterChangeWarningListener warningListener = null;

	/**
	 * Construct attenuation relationship. Initialize parameters and parameter
	 * lists.
	 */
	public AkB_2010_AttenRel(
			final ParameterChangeWarningListener warningListener) {

		// creates exceedProbParam
		super();

		this.warningListener = warningListener;

		initSupportedIntensityMeasureParams();

		// Create an Hash map that links the period with its index
		indexFromPerHashMap = new HashMap<Double, Integer>();
		for (int i = 2; i < AkB2010Constants.PERIOD.length; i++) {
			indexFromPerHashMap.put(new Double(AkB2010Constants.PERIOD[i]),
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
	 * Creates the three supported IM parameters (PGA, SA, and PGV), as well as
	 * the independenParameters of SA (periodParam and dampingParam) and adds
	 * them to the supportedIMParams list. Makes the parameters non-editable.
	 */
	protected final void initSupportedIntensityMeasureParams() {

		// set supported periods for spectral acceleration
		DoubleDiscreteConstraint periodConstraint = new DoubleDiscreteConstraint();
		for (int i = 2; i < AkB2010Constants.PERIOD.length; i++) {
			periodConstraint.addDouble(new Double(AkB2010Constants.PERIOD[i]));
		}
		periodConstraint.setNonEditable();
		// set period param (default is 1s, which is provided by AkB2010 GMPE)
		saPeriodParam = new PeriodParam(periodConstraint);

		// set damping parameter. Empty constructor set damping
		// factor to 5 % (which is the one provided by AkB2010 GMPE)
		saDampingParam = new DampingParam();

		// initialize spectral acceleration parameter (units: g)
		saParam = new SA_Param(saPeriodParam, saDampingParam);
		saParam.setNonEditable();

		// initialize peak ground acceleration parameter (units: g)
		pgaParam = new PGA_Param();
		pgaParam.setNonEditable();

		// initialize peak ground velocity parameter (units: cm/sec)
		pgvParam = new PGV_Param();
		pgvParam.setNonEditable();

		// add the warning listeners
		saParam.addParameterChangeWarningListener(warningListener);
		pgaParam.addParameterChangeWarningListener(warningListener);
		pgvParam.addParameterChangeWarningListener(warningListener);

		// put parameters in the supportedIMParams list
		supportedIMParams.clear();
		supportedIMParams.addParameter(saParam);
		supportedIMParams.addParameter(pgaParam);
		supportedIMParams.addParameter(pgvParam);

	}

	/**
	 * Initialize earthquake rupture parameter (moment magnitude, rake) and add
	 * to eqkRuptureParams list. Makes the parameters non-editable.
	 */
	protected final void initEqkRuptureParams() {

		// moment magnitude
		magParam = new MagParam(AkB2010Constants.MAG_WARN_MIN,
				AkB2010Constants.MAG_WARN_MAX);
		// rake angle (default 0.0 -> strike-slip)
		rakeParam = new RakeParam();

		eqkRuptureParams.clear();
		eqkRuptureParams.addParameter(magParam);
		eqkRuptureParams.addParameter(rakeParam);
	}

	/**
	 * Initialize site parameters (vs30) and adds it to the siteParams list.
	 * Makes the parameters non-editable.
	 */
	protected final void initSiteParams() {

		// vs30 parameters (constrains are not set, default value to 760 m/s)
		vs30Param = new Vs30_Param();

		siteParams.clear();
		siteParams.addParameter(vs30Param);
	}

	/**
	 * Initialize Propagation Effect parameters (Joyner-Boore distance) and adds
	 * them to the propagationEffectParams list. Makes the parameters
	 * non-editable.
	 */
	protected final void initPropagationEffectParams() {

		distanceJBParam = new DistanceJBParameter(
				AkB2010Constants.DISTANCE_JB_WARN_MIN);
		distanceJBParam.addParameterChangeWarningListener(warningListener);
		DoubleConstraint warn = new DoubleConstraint(
				AkB2010Constants.DISTANCE_JB_WARN_MIN,
				AkB2010Constants.DISTANCE_JB_WARN_MAX);
		warn.setNonEditable();
		distanceJBParam.setWarningConstraint(warn);
		distanceJBParam.setNonEditable();

		propagationEffectParams.addParameter(distanceJBParam);
	}

	/**
	 * Initialize other Parameters (standard deviation type, component, sigma
	 * truncation type, sigma truncation level).
	 */
	protected final void initOtherParams() {

		sigmaTruncTypeParam = new SigmaTruncTypeParam();
		sigmaTruncLevelParam = new SigmaTruncLevelParam();

		// stdDevType Parameter
		StringConstraint stdDevTypeConstraint = new StringConstraint();
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_TOTAL);
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_NONE);
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_INTER);
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_INTRA);
		stdDevTypeConstraint.setNonEditable();
		stdDevTypeParam = new StdDevTypeParam(stdDevTypeConstraint);

		// the Component Parameter
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
	protected final void initIndependentParamLists() {

		// params that the mean depends upon
		meanIndependentParams.clear();
		meanIndependentParams.addParameter(magParam);
		meanIndependentParams.addParameter(rakeParam);
		meanIndependentParams.addParameter(distanceJBParam);
		meanIndependentParams.addParameter(vs30Param);

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
	 * Adds the parameter change listeners. This allows to listen to when-ever
	 * the parameter is changed.
	 */
	protected void initParameterEventListeners() {

		magParam.addParameterChangeListener(this);
		rakeParam.addParameterChangeListener(this);
		distanceJBParam.addParameterChangeListener(this);
		vs30Param.addParameterChangeListener(this);
		componentParam.addParameterChangeListener(this);
		stdDevTypeParam.addParameterChangeListener(this);
		saPeriodParam.addParameterChangeListener(this);

	}

	/**
	 * This listens for parameter changes and updates the primitive parameters
	 * accordingly
	 */
	public final void parameterChange(final ParameterChangeEvent e) {

		String pName = e.getParameterName();
		Object val = e.getNewValue();

		if (pName.equals(MagParam.NAME)) {
			mag = ((Double) val).doubleValue();
		} else if (pName.equals(RakeParam.NAME)) {
			rake = ((Double) val).doubleValue();
		} else if (pName.equals(DistanceJBParameter.NAME)) {
			rJB = ((Double) val).doubleValue();
		} else if (pName.equals(Vs30_Param.NAME)) {
			vs30 = ((Double) val).doubleValue();
		} else if (pName.equals(StdDevTypeParam.NAME)) {
			stdDevType = (String) val;
		}
	}

	/**
	 * Allows to reset the change listeners on the parameters
	 */
	public void resetParameterEventListeners() {

		magParam.removeParameterChangeListener(this);
		rakeParam.removeParameterChangeListener(this);
		distanceJBParam.removeParameterChangeListener(this);
		vs30Param.removeParameterChangeListener(this);
		stdDevTypeParam.removeParameterChangeListener(this);
		this.initParameterEventListeners();
	}

	/**
	 * This sets the eqkRupture related parameters (moment magnitude, rake
	 * angle) based on the eqkRupture passed in. The internally held eqkRupture
	 * object is also set as that passed in.
	 */
	public final void setEqkRupture(final EqkRupture eqkRupture) {

		magParam.setValueIgnoreWarning(new Double(eqkRupture.getMag()));
		if (!Double.isNaN(eqkRupture.getAveRake())) {
			rakeParam.setValue(eqkRupture.getAveRake());
		}
		this.eqkRupture = eqkRupture;
		setPropagationEffectParams();
	}

	/**
	 * This sets the site-related parameter (vs30) based on what is in the Site
	 * object passed in. This also sets the internally held Site object as that
	 * passed in.
	 */
	public final void setSite(final Site site) {

		vs30Param.setValueIgnoreWarning((Double) site.getParameter(
				Vs30_Param.NAME).getValue());
		this.site = site;
		setPropagationEffectParams();
	}

	/**
	 * This calculates the JB Distance propagation effect parameter based on the
	 * current site and eqkRupture.
	 */
	protected void setPropagationEffectParams() {

		if ((this.site != null) && (this.eqkRupture != null)) {
			distanceJBParam.setValue(eqkRupture, site);
		}
	}

	/**
	 * Set period index.
	 */
	protected final void setPeriodIndex() {
		if (im.getName().equalsIgnoreCase(PGV_Param.NAME)) {
			iper = 0;
		} else if (im.getName().equalsIgnoreCase(PGA_Param.NAME)) {
			iper = 1;
		} else {
			iper = ((Integer) indexFromPerHashMap.get(saPeriodParam.getValue()))
					.intValue();
		}
	}

	/**
	 * Compute mean.
	 */
	public double getMean() {

		if (rJB > USER_MAX_DISTANCE) {
			return VERY_SMALL_MEAN;
		} else {
			setPeriodIndex();
			return getMean(iper, mag, rJB, vs30, rake);
		}
	}

	/**
	 * Compute standard deviation.
	 */
	public final double getStdDev() {

		setPeriodIndex();
		return getStdDev(iper, stdDevType);
	}

	/**
	 * Allows the user to set the default parameter values for the selected
	 * Attenuation Relationship.
	 */
	public final void setParamDefaults() {

		magParam.setValueAsDefault();
		rakeParam.setValueAsDefault();
		distanceJBParam.setValueAsDefault();
		vs30Param.setValueAsDefault();
		saPeriodParam.setValueAsDefault();
		saDampingParam.setValueAsDefault();
		saParam.setValueAsDefault();
		pgaParam.setValueAsDefault();
		pgvParam.setValueAsDefault();
		stdDevTypeParam.setValueAsDefault();
		sigmaTruncTypeParam.setValueAsDefault();
		sigmaTruncLevelParam.setValueAsDefault();
		componentParam.setValueAsDefault();
	}

	/**
	 * Get the name of this IMR.
	 */
	public final String getName() {
		return NAME;
	}

	/**
	 * Returns the Short Name of each AttenuationRelationship
	 * 
	 */
	public final String getShortName() {
		return SHORT_NAME;
	}

	/**
	 * Compute mean (natural logarithm of median ground motion).
	 */
	public double getMean(int iper, double mag, final double rJB,
			final double vs30, final double rake) {

		double logY = Double.NaN;

		double[] b = setConstants(iper);

		int[] soilTerms = setSoilTerms(vs30);

		int[] faultStyleTerms = setFaultStyleTerms(rake);

		logY = b[0] + b[1] * mag + b[2] * mag * mag + (b[3] + b[4] * mag)
				* Math.log10(Math.sqrt(rJB * rJB + b[5] * b[5])) + b[6]
				* soilTerms[0] + b[7] * soilTerms[1] + b[8]
				* faultStyleTerms[0] + b[9] * faultStyleTerms[1];

		// convert from log10 to natural logarithm
		logY = logY * AkB2010Constants.LOG10_2_LN;

		// convert PGA and SA from cm/s2 to g
		if (iper == 0) {
			logY = Math.log(Math.exp(logY));
		} else {
			logY = Math.log(Math.exp(logY)
					* AkB2010Constants.CMS2_TO_G_CONVERSION_FACTOR);
		}
		
		// apply correction for T = 4
		if (AkB2010Constants.PERIOD[iper] == 4.00){
			logY /= AkB2010Constants.T3sec_TO_T4sec_DECAYFACTOR;
		}

		return logY;
	}

	private int[] setFaultStyleTerms(final double rake) {
		int[] faultStyleTerms = new int[] { 0, 0 };
		boolean normal = rake > -150.0 && rake < -30.0;
		boolean reverse = rake > 30.0 && rake < 150.0;
		if (normal) {
			faultStyleTerms[0] = 1;
		}
		if (reverse) {
			faultStyleTerms[1] = 1;
		}
		return faultStyleTerms;
	}

	private int[] setSoilTerms(final double vs30) {
		int[] soilTerms = new int[] { 0, 0 };

		if (vs30 < AkB2010Constants.SOFT_SOIL_UPPER_BOUND) {
			soilTerms[0] = 1;
		}
		if (vs30 >= AkB2010Constants.SOFT_SOIL_UPPER_BOUND
				&& vs30 <= AkB2010Constants.STIFF_SOIL_UPPER_BOUND) {
			soilTerms[1] = 1;
		}
		return soilTerms;
	}

	private double[] setConstants(int iper) {
		double[] b = new double[10];

		b[0] = AkB2010Constants.b1[iper];

		b[1] = AkB2010Constants.b2[iper];

		b[2] = AkB2010Constants.b3[iper];

		b[3] = AkB2010Constants.b4[iper];

		b[4] = AkB2010Constants.b5[iper];

		b[5] = AkB2010Constants.b6[iper];

		b[6] = AkB2010Constants.b7[iper];

		b[7] = AkB2010Constants.b8[iper];

		b[8] = AkB2010Constants.b9[iper];

		b[9] = AkB2010Constants.b10[iper];
		return b;
	}

	public double getStdDev(int iper, String stdDevType) {
		if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_NONE))
			return 0;
		else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_INTER))
			return AkB2010Constants.LOG10_2_LN
					* AkB2010Constants.INTER_EVENT_STD[iper];
		else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_INTRA))
			return AkB2010Constants.LOG10_2_LN
					* AkB2010Constants.INTRA_EVENT_STD[iper];
		else
			return AkB2010Constants.LOG10_2_LN
					* AkB2010Constants.TOTAL_STD[iper];
	}

	/**
	 * This provides a URL where more info on this model can be obtained.
	 * Currently returns null because no URL has been created
	 */
	public URL getInfoURL() throws MalformedURLException {
		return null;
	}

}
