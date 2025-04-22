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
	enum ebranch { gapping, precap, postcap, residual, failing, failed, peakOriented, pinching, peakPassing };
	SmoothIMK(int tag,
		double pd1, double pf1,
		double pd2, double pf2,
		double pd3, double pf3,
		double pdu, 
		double nd1, double nf1,
		double nd2, double nf2,
		double nd3, double nf3,
		double ndu,
		double gamaS, double cS,
		double gamaUnloadE, double cUnloadE,
		double r0, double r1, double r2,
		int cyclicRule,
		double pinchXPos, double pinchYPos, double sigPenetFacPos,
		double pinchXNeg, double pinchYNeg, double sigPenetFacNeg,
		double sigInit, double gap);

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

	double getInitYieldStrain() { return pd1_; }
	virtual void resetEnergy(void) { EnergyP = 0; }
	virtual Response* setResponse(const char** argv, int argc,
		OPS_Stream& theOutputStream);
	virtual int getResponse(int responseID, Information& matInformation);

protected:

private:
	double pf1, pd1, pf2, pd2, pf3, pd3, pdu;
	double pd1_, pd2_, pd3_, pdu_, nd1_, nd2_, nd3_, ndu_;
	double nf1, nd1, nf2, nd2, nf3, nd3, ndu;
	double r0, r1, r2;
	double gap;
	double sigPenetFacP, sigPenetFacN;
	int cyclicRule;  // 1:bilinear, 2:pinched, 3:peak-oriented
	double FailEnergS, FailEnergUnloadE, cS, cUnloadE;			//damage parameters
	double pinchXPos, pinchXNeg, pinchYPos, pinchYNeg;
	double FydP, FydN;		//Pos and Neg Fy's affected by damage
	double ExcurEnergy;
	double EnergyP; //by SAJalali
	double sigini; // initial 
	
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
	ebranch branch, branchP;
	// hstv : STEEL HISTORY VARIABLES   
	double epsP;  //  = strain at previous converged step
	double sigP;  //  = stress at previous converged step
	double eP;    //   stiffness modulus at last converged step;
	double slopeRat, slopeRatP;
	double epsmin;
	double epsmax;
	double epsLimit;
	double epss0;
	double sigs0;
	double epsr;
	double sigr;
	double sig;
	double e;
	double eps;   //  = strain at current step
	double EshP, EshN, E0p, E0n, EunloadP, EunloadN;
	double R0, R0P;
	double FcP, FcN;
	double FrP, FrN;
	bool onEnvelope, onEnvelopeP;
	bool initiated, initiatedP;
	double epsPl, epsPlP;
	void updateDamage();
	void updateAsymptote();
	void bilinAsymptote();
	void peakOrientedAsymptote();
	void pinchedAsymptote();
	void changeBranch(bool isReturning);
	void gappedBilinAsymptote();
	ebranch nextBranch(ebranch branch);
	void getEnvelope(double eps, double& targStress, SmoothIMK::ebranch& targBranch, double& k, double& limitEps);
	void computeR0(double k1, double k2, double E1, double dy);
};


#endif

