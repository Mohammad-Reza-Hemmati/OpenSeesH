/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.0 $                                                                     
// Written: 
// IMK material with exponential transition curves and multiple cyclic modes (Bilinear, poinched and peak-oriented)


#ifndef SmoothIMK_h
#define SmoothIMK_h

#include <UniaxialMaterial.h>
class SmoothIMK : public UniaxialMaterial
{
public:
	//enum ebranch { gapping, precap, postcap, residual, failing, failed, peakOriented, pinching, peakPassing };
	SmoothIMK(int tag,
		std::vector<double> pd, std::vector<double> pf,
		std::vector<double> nd, std::vector<double> nf,
		double gamaS, double cS,
		double gamaUnloadE, double cUnloadE,
		double r0, double r1, double r2,
		int cyclicRuleP, double gapP, double alphaPinchPos, double pinchYPos,
		double betaPinchPos, double epsPinchPos, double sigPenetFacPos, double bilinEndAmpPos, double unloadingStiffFacPos,
		int cyclicRuleN, double gapN, double alphaPinchNeg, double pinchYNeg,
		double betaPinchNeg, double epsPinchNeg, double sigPenetFacNeg, double bilinEndAmpNeg, double unloadingStiffFacNeg,
		double sigInit);

	SmoothIMK(void);
	virtual ~SmoothIMK();


	const char* getClassType(void) const { return "SmoothIMK"; };

	double getInitialTangent(void);
	UniaxialMaterial* getCopy(void);

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);

	virtual double getEnergy() { return EnergyP; };

	double getInitYieldStrain() { return pd[0]; }
	virtual void resetEnergy(void) { EnergyP = 0; }
	virtual Response* setResponse(const char** argv, int argc,
		OPS_Stream * theOutputStream);
	virtual int getResponse(int responseID, Information& matInformation);

protected:

private:
	std::vector<double> pd, pf, nd, nf, gpd, gnd; //gpd: gapped pd
	double r0, r1, r2;
	double gapP, gapN;
	double bilinEndAmpP, bilinEndAmpN;	//strain amplitude limit at which bilin behavior ends and either pinched or peak-oriented rules start based on the cyclicRule
	double sigPenetFacP, sigPenetFacN;
	int cyclicRuleP, cyclicRuleN;  // 1:bilinear, 2:pinched, 3:peak-oriented
	double FailEnergS, FailEnergUnloadE, cS, cUnloadE;			//damage parameters
	double alphaPinchPos, alphaPinchNeg, betaPinchPos, epsPinchPos, pinchYPos, pinchYNeg, betaPinchNeg, epsPinchNeg;
	//double FydP, FydN;		//Pos and Neg Fy's affected by damage
	double ExcurEnergy;
	double EnergyP; //by SAJalali
	double sigini; // initial 
	double unloadingStiffFacPos, unloadingStiffFacNeg;
	//HISTORY VARIABLES
	double epsminP; //  = hstvP(1) : max eps in compression
	double epsmaxP; //  = hstvP(2) : max eps in tension
	double epsLimitP;  //  = hstvP(3) : plastic excursion
	double epss0P;  //  = hstvP(4) : eps at asymptotes intersection
	double sigs0P;  //  = hstvP(5) : sig at asymptotes intersection
	double epssrP;  //  = hstvP(6) : eps at last inversion point
	double sigsrP;  //  = hstvP(7) : sig at last inversion point
	//int    stat, statP;    //  = hstvP(8) : index for loading/unloading
	bool isPosDir, isPosDirP;
	int branch, branchP; 
	// 0: initial on-envelope branch,
	// 1~n: segment-wise on-envelope,
	// 10011: initial gapped loading,
	// 10012: gapped unloading,
	// 1002: pinching,
	// 1003: peak-oriented,
	// 1005: peak-passing,
	// 1010: failed

	double eps, epsP;  //  = strain at trial and last converged steps
	double sig, sigP;  //  = stress at traial and converged steps
	double e, eP;    //   stiffness modulus at trial and last converged stepa;
	double slopeRat, slopeRatP;
	double epsPl, epsPlP;	//cumulative plastic strain at trial and last converged steps
	double epsmin;
	double epsmax;
	double epsLimit;
	double epss0;
	double sigs0;
	double epsr;
	double sigr;
	double E0p, E0n, EunloadP, EunloadN;
	//double EshP, EshN;
	double R0, R0P;
	std::vector<double> fpDmgd, fnDmgd;
	int FyIndP,FyIndN, FcIndP, FcIndN;
	int FrIndP, FrIndN;
	bool onEnvelope, onEnvelopeP;
	bool initiated, initiatedP;
	void updateDamage();
	void updateAsymptote();
	void changeBranch(bool isReturning);
	int nextBranch(int branch, bool shouldPinch);
	void getEnvelope(double eps, bool pSide, double& targStress, int& targBranch, double& k, double& limitEps);
	void computeR0(double k1, double k2, double E1, double dy);
};


#endif

