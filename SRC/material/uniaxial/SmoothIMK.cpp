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

// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SmoothIMK.cpp,v $

// Written: S. A. Jalali 10/2019
// Adding Cyclic and in-cycle deterioration modes to steel02 UniaxialMaterial

#include <math.h>

#include <stdlib.h>
#include <SmoothIMK.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include "classTags.h"

static int numThisCall = 0;
void printSyntax()
{
	opserr << "------ SmoothIMK unaxialMaterial -------\n";
	opserr << "-------Syntax:\n";
	opserr << "-------UniaxialMaterial SmoothIMK $matTag\n";
	opserr << "									 pd1, pf1, pd2, pf2, pd3, pf3, pdu\n";
	opserr << "                  <nd1, nf1, nd2, nf2, nd3, nf3, ndu>\n";
	opserr << "                  <-sigInit sigInit>\n";
	opserr << "                  <-gap gap>\n";
	opserr << "                  <-deterioration gamaS cS gamaK cK>\n";
	opserr << "                  <-transition r0 <r1 r2>> \n";
	opserr << "                  <-peakOriented <stressPenetFacPos <stressPenetFacNeg>>> \n";
	opserr << "                  <-pinched <stressPenetFacPos pinchXPos pinchYPos <stressPenetFacNeg pinchXNeg pinchYNeg>>> \n";
}
void*
OPS_SmoothIMK()
{
	if (numThisCall == 0) {
		numThisCall = 1;
	}
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	int numData = 1;
	int tag = 0;
	if (OPS_GetIntInput(&numData, &tag) != 0) {
		opserr << "WARNING invalid uniaxialMaterial SmoothIMK tag" << endln;
		printSyntax();
		return 0;
	}

	if (OPS_GetNumRemainingInputArgs() < 7) {
		opserr << "Invalid SmoothIMK #args for: " << tag << " see the syntax" << endln;
		printSyntax();
		return 0;
	}
	double pf1, pd1, pf2, pd2, pf3, pd3, pdu;
	double nf1, nd1, nf2, nd2, nf3, nd3, ndu;
	//default parameters:
	double gamaS = 1e6, cS = 1, gamaUE = 1e6, cUE = 1, r0 = 1.4, r1 = 26.6, r2 = 0.84, sigInit = 0;
	double pinchXPos = 0.5, pinchXNeg = 0.5, pinchYPos = 0.2, pinchYNeg = 0.2;
	double sigPenetFacP = 0, sigPenetFacN = 0;
	double gap = 0;
	int rule = 1;

	if (OPS_GetDoubleInput(&numData, &pd1) != 0) {
		opserr << "SmoothIMK:: invalid pd1 for material : " << tag << endln;
		printSyntax();
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pf1) != 0) {
		opserr << "SmoothIMK:: invalid pf1 for material : " << tag << endln;
		printSyntax();
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pd2) != 0) {
		opserr << "SmoothIMK:: invalid pd2 for material : " << tag << endln;
		printSyntax();
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pf2) != 0) {
		opserr << "SmoothIMK:: invalid pf2 for material : " << tag << endln;
		printSyntax();
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pd3) != 0) {
		opserr << "SmoothIMK:: invalid pd3 for material : " << tag << endln;
		printSyntax();
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pf3) != 0) {
		opserr << "SmoothIMK:: invalid pf3 for material : " << tag << endln;
		printSyntax();
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &pdu) != 0) {
		opserr << "SmoothIMK:: invalid pdu for material : " << tag << endln;
		printSyntax();
		return 0;
	}
	if (pd1 >= pd2 || pd2 >= pd3 || pd3 > pdu)
	{
		opserr << "SmoothIMK:: positive strain limits should meet: pd1 < pd2 < pd3 <= pdu; for material : " << tag << endln;
		printSyntax();
		return 0;
	}
	if (OPS_GetDoubleInput(&numData, &nd1) != 0) {
		//OPS_ResetCurrentInputArg(-1);
		nd1 = -pd1, nf1 = -pf1, nf2 = -pf2, nd2 = -pd2, nf3 = -pf3, nd3 = -pd3, ndu = -pdu;
	}
	else {
		if (OPS_GetDoubleInput(&numData, &nf1) != 0) {
			opserr << "SmoothIMK:: invalid nf1 for material : " << tag << endln;
			printSyntax();
			return 0;
		}
		if (OPS_GetDoubleInput(&numData, &nd2) != 0) {
			opserr << "SmoothIMK:: invalid nd2 for material : " << tag << endln;
			printSyntax();
			return 0;
		}
		if (OPS_GetDoubleInput(&numData, &nf2) != 0) {
			opserr << "SmoothIMK:: invalid nf2 for material : " << tag << endln;
			printSyntax();
			return 0;
		}
		if (OPS_GetDoubleInput(&numData, &nd3) != 0) {
			opserr << "SmoothIMK:: invalid nd3 for material : " << tag << endln;
			printSyntax();
			return 0;
		}
		if (OPS_GetDoubleInput(&numData, &nf3) != 0) {
			opserr << "SmoothIMK:: invalid nf3 for material : " << tag << endln;
			printSyntax();
			return 0;
		}
		if (OPS_GetDoubleInput(&numData, &ndu) != 0) {
			opserr << "SmoothIMK:: invalid ndu for material : " << tag << endln;
			printSyntax();
			return 0;
		}
		nd1 = -fabs(nd1); // Ensure nf1 is negative
		nf1 = -fabs(nf1); // Ensure nf1 is negative
		nd2 = -fabs(nd2); // Ensure nd2 is negative
		nf2 = -fabs(nf2); // Ensure nf2 is negative
		nd3 = -fabs(nd3); // Ensure nd3 is negative
		nf3 = -fabs(nf3); // Ensure nf3 is negative
		ndu = -fabs(ndu); // Ensure ndu is negative
		if (nd1 <= nd2 || nd2 <= nd3 || nd3 < ndu)
		{
			opserr << "SmoothIMK:: negative strain limits should meet: nd1 > nd2 > nd3 >= ndu; for material : " << tag << endln;
			printSyntax();
			return 0;
		}
	}
	while (OPS_GetNumRemainingInputArgs() > 0)
	{
		const char* option = OPS_GetString();
		if (strcmp(option, "-sigInit") == 0) {
			if (OPS_GetDoubleInput(&numData, &sigInit) != 0) {
				opserr << "SmoothIMK:: invalid sigInit for material : " << tag << endln;
				return 0;
			}
		}
		else if (strcmp(option, "-gap") == 0) {
			if (OPS_GetDoubleInput(&numData, &gap) != 0) {
				opserr << "SmoothIMK:: invalid gap for material : " << tag << endln;
				return 0;
			}
		}
		else if (strcmp(option, "-deterioration") == 0) {
			if (OPS_GetDoubleInput(&numData, &gamaS) != 0) {
				opserr << "SmoothIMK:: invalid gamaS for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &cS) != 0) {
				opserr << "SmoothIMK:: invalid cS for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &gamaUE) != 0) {
				opserr << "SmoothIMK:: invalid gamaUloadE for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &cUE) != 0) {
				opserr << "SmoothIMK:: invalid cUloadE for material : " << tag << endln;
				return 0;
			}
		}
		else if (strcmp(option, "-transition") == 0) {
			if (OPS_GetDoubleInput(&numData, &r0) != 0) {
				opserr << "SmoothIMK:: invalid r0 for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &r1) != 0) {
				//OPS_ResetCurrentInputArg(-1);
				r1 = 0, r2 = 0;
			}
			else {
				if (OPS_GetDoubleInput(&numData, &r2) != 0) {
					opserr << "SmoothIMK:: invalid r2 for material : " << tag << endln;
					return 0;
				}
			}
		}
		else if (strcmp(option, "-peakOriented") == 0) {
			rule = 2;
			if (OPS_GetDoubleInput(&numData, &sigPenetFacP) != 0)
				//OPS_ResetCurrentInputArg(-1);
				;
			else
				if (OPS_GetDoubleInput(&numData, &sigPenetFacN) != 0)
				{
					//OPS_ResetCurrentInputArg(-1);
					sigPenetFacN = sigPenetFacP;
				}
		}
		else if (strcmp(option, "-pinched") == 0) {
			rule = 3;
			if (OPS_GetDoubleInput(&numData, &sigPenetFacP) != 0) {
				opserr << "SmoothIMK:: invalid sigPenetFacPos for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &pinchXPos) != 0) {
				opserr << "SmoothIMK:: invalid pinchXPos for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &pinchYPos) != 0) {
				opserr << "SmoothIMK:: invalid pinchYPos for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &sigPenetFacN) != 0) {
				//OPS_ResetCurrentInputArg(-1);
				sigPenetFacN = sigPenetFacP;
				pinchXNeg = pinchXPos;
				pinchYNeg = pinchYPos;
			}
			else {
				if (OPS_GetDoubleInput(&numData, &pinchXNeg) != 0) {
					opserr << "SmoothIMK:: invalid pinchXNeg for material : " << tag << endln;
					return 0;
				}
				if (OPS_GetDoubleInput(&numData, &pinchYNeg) != 0) {
					opserr << "SmoothIMK:: invalid pinchYNeg for material : " << tag << endln;
					return 0;
				}
			}
		}
		else {
			if (strcmp(option, "Invalid String Input!") == 0)
			{
				//OPS_ResetCurrentInputArg(-1);
				char data[128];
				OPS_GetStringFromAll(data, 128);
				opserr << "SmoothIMK::Error: expected string switch but got non-string input: " << data << " for material with tag: " << tag << endln;
			}
			else
				opserr << "SmoothIMK::Error: Invalid option: " << option << endln;
			return 0;
		}
	}
	// Parsing was successful, allocate the material
	theMaterial = new SmoothIMK(tag, pd1, pf1, pd2, pf2, pd3, pf3, pdu, nd1, nf1, nd2, nf2, nd3, nf3, ndu,
		gamaS, cS, gamaUE, cUE, r0, r1, r2, rule, pinchXPos, pinchYPos, sigPenetFacP,
		pinchXNeg, pinchYNeg, sigPenetFacN, sigInit, gap);


	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type SmoothIMK\n";
		return 0;
	}

	return theMaterial;
}

SmoothIMK::SmoothIMK(int tag,
	double _pd1, double _pf1,
	double _pd2, double _pf2,
	double _pd3, double _pf3,
	double _pdu,
	double _nd1, double _nf1,
	double _nd2, double _nf2,
	double _nd3, double _nf3,
	double _ndu,
	double _gamaS, double _cS,
	double _gamaUE, double _cUE,
	double _r0, double _r1, double _r2,
	int _cyclicRule,
	double _pinchXPos, double _pinchYPos, double _sigPenetFacP,
	double _pinchXNeg, double _pinchYNeg, double _sigPenetFacN,
	double sigInit, double _gap) :
	UniaxialMaterial(tag, MAT_TAG_SmoothIMK),
	pd1_(_pd1), pf1(_pf1),
	pd2_(_pd2), pf2(_pf2),
	pd3_(_pd3), pf3(_pf3),
	nd1_(_nd1), nf1(_nf1),
	nd2_(_nd2), nf2(_nf2),
	nd3_(_nd3), nf3(_nf3),
	pdu_(_pdu), ndu_(_ndu),
	cS(_cS), FailEnergS(_gamaS* pf1* pd1_),
	cUnloadE(_cUE), FailEnergUnloadE(_gamaUE* pf1* pd1_),
	r0(_r0), r1(_r1), r2(_r2),
	cyclicRule(_cyclicRule),
	pinchXPos(_pinchXPos), pinchYPos(_pinchYPos),
	pinchXNeg(_pinchXNeg), pinchYNeg(_pinchYNeg),
	sigPenetFacP(_sigPenetFacP), sigPenetFacN(_sigPenetFacN),
	sigini(sigInit), gap(_gap)
{
	revertToStart();
}

SmoothIMK::SmoothIMK(void) :
	UniaxialMaterial(0, MAT_TAG_SmoothIMK),
	pd1_(0), pf1(0),
	pd2_(0), pf2(0),
	pd3_(0), pf3(0),
	nd1_(0), nf1(0),
	nd2_(0), nf2(0),
	nd3_(0), nf3(0),
	pdu_(0), ndu_(0),
	cS(0), FailEnergS(0),
	cUnloadE(0), FailEnergUnloadE(0),
	r0(0), r1(0), r2(0),
	cyclicRule(1),
	sigini(0)
{
	//check inputs
	//sigIni and gap should not be given simultaneously
	//pd1 > 0, pdu >= pd3 > pd2 > pd1
	//nd1 < 0, ndu <= nd3 < nd2 < pd1
	revertToStart();
}

SmoothIMK::~SmoothIMK(void)
{
	// Does nothing
}

UniaxialMaterial*
SmoothIMK::getCopy(void)
{
	SmoothIMK* theCopy = new SmoothIMK(this->getTag(), pd1_, pf1, pd2_, pf2, pd3_,
		pf3, pdu_, nd1_, nf1, nd2_, nf2, nd3_, nf3, ndu_, FailEnergS / pf1 / pd1_, cS, FailEnergUnloadE / pf1 / pd1_, cUnloadE, r0, r1, r2, cyclicRule,
		pinchXPos, pinchYPos, sigPenetFacP, pinchXNeg, pinchYNeg, sigPenetFacN, sigini, gap);
	theCopy->revertToStart();
	return theCopy;
}

double
SmoothIMK::getInitialTangent(void)
{
	return E0p;
}

double SmoothIMK::getStrain(void)
{
	return eps;
}

double SmoothIMK::getStress(void)
{
	return sig;
}

double SmoothIMK::getTangent(void)
{
	return e;
}

int
SmoothIMK::setParameter(const char** argv, int argc, Parameter& param)
{
	return -1;
}

int
SmoothIMK::updateParameter(int parameterID, Information& info)
{
	return -1;
}

int
SmoothIMK::commitState(void)
{
	epsminP = epsmin;
	epsmaxP = epsmax;
	epsLimitP = epsLimit;
	epss0P = epss0;
	sigs0P = sigs0;
	epssrP = epsr;
	sigsrP = sigr;
	slopeRatP = slopeRat;
	onEnvelopeP = onEnvelope;
	epsPlP = epsPl;
	updateDamage();
	isPosDirP = isPosDir;
	branchP = branch;
	eP = e;
	sigP = sig;
	epsP = eps;
	R0P = R0;
	initiatedP = initiated;
	return 0;
}

int
SmoothIMK::revertToLastCommit(void)
{
	epsmin = epsminP;
	epsmax = epsmaxP;
	epsLimit = epsLimitP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	branch = branchP;
	isPosDir = isPosDirP;
	slopeRat = slopeRatP;
	R0 = R0P;
	e = eP;
	sig = sigP;
	eps = epsP;
	onEnvelope = onEnvelopeP;
	epsPl = epsPlP;
	initiated = initiatedP;
	return 0;
}

int
SmoothIMK::revertToStart(void)
{
	EnergyP = 0;	//by SAJalali
	pd1 = pd1_;
	pd2 = pd2_;
	pd3 = pd3_;
	pdu = pdu_;
	nd1 = nd1_;
	nd2 = nd2_;
	nd3 = nd3_;
	ndu = ndu_;
	E0p = pf1 / pd1_;
	E0n = nf1 / nd1_;
	e = eP = E0p;
	onEnvelope = onEnvelopeP = true;
	branchP = branch = precap;
	FydP = pf1;
	FydN = nf1;
	FcP = pf2;
	FcN = nf2;
	EshP = (FcP - FydP) / (pd2 - pd1);
	EshN = (FcN - FydN) / (nd2 - nd1);
	FrP = pf3;
	FrN = nf3;
	epsP = sigini / E0p;
	sigP = sigini;
	sig = 0.0;
	eps = 0.0;
	EunloadP = E0p;
	EunloadN = E0n;
	isPosDirP = isPosDir = true;
	epsmaxP = pd1;
	epsminP = nd1;
	epsLimitP = pd1;
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;
	R0P = r0;
	epsPl = epsPlP = 0;
	ExcurEnergy = 0;
	slopeRat = slopeRatP = 0;
	initiated = initiatedP = false;
	return 0;
}

int
SmoothIMK::sendSelf(int commitTag, Channel& theChannel)
{
	static Vector data(60);
	int n = -1;
	data(n++) = this->getTag();	//0
	data(n++) = pd1_;
	data(n++) = pd2_;
	data(n++) = pd3_;
	data(n++) = pdu_;
	data(n++) = nd1_;
	data(n++) = nd2_;
	data(n++) = nd3_;
	data(n++) = ndu_;
	data(n++) = pf1;
	data(n++) = pf2;
	data(n++) = pf3;
	data(n++) = nf1;
	data(n++) = nf2;
	data(n++) = nf3;
	data(n++) = pinchXPos;
	data(n++) = pinchYPos;
	data(n++) = sigPenetFacP;
	data(n++) = pinchXNeg;
	data(n++) = pinchYNeg;
	data(n++) = sigPenetFacN;
	data(n++) = FailEnergS;
	data(n++) = cS;
	data(n++) = FailEnergUnloadE;
	data(n++) = cUnloadE;
	data(n++) = sigini;
	data(n++) = epsP;
	data(n++) = sigP;
	data(n++) = eP;
	data(n++) = EnergyP;
	data(n++) = epsminP;
	data(n++) = epsmaxP;
	data(n++) = epsLimitP;
	data(n++) = epss0P;
	data(n++) = sigs0P;
	data(n++) = epssrP;
	data(n++) = sigsrP;
	data(n++) = isPosDirP;
	data(n++) = branchP;
	data(n++) = FydP;
	data(n++) = FydN;
	data(n++) = FcP;
	data(n++) = FrP;
	data(n++) = FcN;
	data(n++) = FrN;
	data(n++) = epsPlP;
	data(n++) = ExcurEnergy;
	data(n++) = slopeRatP;
	data(n++) = onEnvelopeP;
	data(n++) = initiatedP;
	data(n++) = R0P;
	data(n++) = pd1;
	data(n++) = pd2;
	data(n++) = pd3;
	data(n++) = pdu;
	data(n++) = nd1;
	data(n++) = nd2;
	data(n++) = nd3;
	data(n++) = ndu;
	data(n++) = gap;

	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SmoothIMK::sendSelf() - failed to sendSelf\n";
		return -1;
	}
	return 0;
}

int
SmoothIMK::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	static Vector data(60);	//editted by SAJalali for energy

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SmoothIMK::recvSelf() - failed to recvSelf\n";
		return -1;
	}
	int n = -1;
	this->setTag(data(n++));
	pd1_ = data(n++);
	pd2_ = data(n++);
	pd3_ = data(n++);
	pdu_ = data(n++);
	nd1_ = data(n++);
	nd2_ = data(n++);
	nd3_ = data(n++);
	ndu_ = data(n++);
	pf1 = data(n++);
	pf2 = data(n++);
	pf3 = data(n++);
	nf1 = data(n++);
	nf2 = data(n++);
	nf3 = data(n++);
	pinchXPos = data(n++);
	pinchYPos = data(n++);
	sigPenetFacP = data(n++);
	pinchXNeg = data(n++);
	pinchYNeg = data(n++);
	sigPenetFacN = data(n++);
	FailEnergS = data(n++);
	cS = data(n++);
	FailEnergUnloadE = data(n++);
	cUnloadE = data(n++);
	sigini = data(n++);
	epsP = data(n++);
	sigP = data(n++);
	eP = data(n++);
	EnergyP = data(n++);
	epsminP = data(n++);
	epsmaxP = data(n++);
	epsLimitP = data(n++);
	epss0P = data(n++);
	sigs0P = data(n++);
	epssrP = data(n++);
	sigsrP = data(n++);
	isPosDirP = data(n++);
	branchP = (ebranch)data(n++);
	FydP = data(n++);
	FydN = data(n++);
	FcP = data(n++);
	FrP = data(n++);
	FcN = data(n++);
	FrN = data(n++);
	epsPlP = data(n++);
	ExcurEnergy = data(n++);
	slopeRatP = data(n++);
	onEnvelopeP = data(n++);
	initiatedP = data(n++);
	R0P = data(n++); //49
	pd1 = data(n++);
	pd2 = data(n++);
	pd3 = data(n++);
	pdu = data(n++);
	nd1 = data(n++);
	nd2 = data(n++);
	nd3 = data(n++);
	ndu = data(n++);
	gap = data(n++);
	e = eP;
	sig = sigP;
	eps = epsP;

	return 0;
}

void
SmoothIMK::Print(OPS_Stream& s, int flag)
{
	//    s << "SmoothIMK:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
	const char* endStr = endln;
	if (flag == OPS_PRINT_PRINTMODEL_JSON)
	{
		endStr = "";
		s << "\t\t\t{";
	}
	s << "\"name \": \"" << this->getTag() << "\", " << endStr;
	s << "\"type\": \"SmoothIMK\", " << endStr;
	s << "\"pd1 \": " << pd1 << ", " << endStr;
	s << "\"pd2 \": " << pd2 << ", " << endStr;
	s << "\"pd3 \": " << pd3 << ", " << endStr;
	s << "\"nd1 \": " << nd1 << ", " << endStr;
	s << "\"nd2 \": " << nd2 << ", " << endStr;
	s << "\"nd3 \": " << nd3 << ", " << endStr;
	s << "\"pf1 \": " << pf1 << ", " << endStr;
	s << "\"pf2 \": " << pf2 << ", " << endStr;
	s << "\"pf3 \": " << pf3 << ", " << endStr;
	s << "\"nf1 \": " << nf1 << ", " << endStr;
	s << "\"nf2 \": " << nf2 << ", " << endStr;
	s << "\"nf3 \": " << nf3 << ", " << endStr;
	s << "\"sigini \": " << sigini << endStr;
	if (flag == OPS_PRINT_PRINTMODEL_JSON)
		s << "}";

}


Response* SmoothIMK::setResponse(const char** argv, int argc, OPS_Stream& theOutput)
{
	Response* theResponse = UniaxialMaterial::setResponse(argv, argc, theOutput);
	if (theResponse != 0)
		return theResponse;
	if (strcmp(*argv, "branch") == 0)
	{
		theResponse = new MaterialResponse(this, 11, 0);
		return theResponse;
	}
	return 0;
}

int SmoothIMK::getResponse(int responseID, Information& matInformation)
{
	int res = 0;
	res = UniaxialMaterial::getResponse(responseID, matInformation);
	if (res == 0)
		return 0;
	switch (responseID)
	{
	case 11:
		matInformation.setInt(branchP);
		return 0;
	default:
		return -1;
	}
}


void SmoothIMK::updateDamage()
{
	if (sigP * sig < 0)
	{
		if (sigP > 0)
		{
			double zeroSigEps = epsP - sigP / E0p;
			double dE = 0.5 * sigP * (zeroSigEps - epsP);
			EnergyP += dE;
			if (EnergyP < 0) EnergyP = 0.;
			ExcurEnergy += dE;
			if (ExcurEnergy < 0) ExcurEnergy = 0.;
			if (branch == failing)
			{
				FydP = FrP;
				FcP = FrP;
				EshP = (FcP - FydP) / (pd2 - pd1);
			}
			else if (branch == failed)
			{
				FydP = 0;
				FcP = FrP;
				EshP = (FcP - FydP) / (pd2 - pd1);
			}
			else
			{
				double beta = 0;
				if (EnergyP > FailEnergS)
					beta = 1;
				else
					beta = pow(ExcurEnergy / (FailEnergS - EnergyP), cS);
				if (beta > 0.9999)
				{
					opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Complete Strength loss\n" << endln;
					beta = 0.9999;
				}
				FydP = (1. - beta) * FydP + beta * pf3 / pf1 * FydP;
				FcP = pf2 - pf1 + FydP;
				EshP = (FcP - FydP) / (pd2 - pd1);
				//FrP = FcP - pf2 + pf3;
				//if (FrP < 0)
				//	FrP = 0;

				if (EnergyP > FailEnergUnloadE)
					beta = 1;
				else
					beta = pow(ExcurEnergy / (FailEnergUnloadE - EnergyP), cUnloadE);
				if (beta > 0.9999)
				{
					opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Complete Unloading Stiffness loss\n" << endln;
					beta = 0.9999;
				}
				EunloadP = (1. - beta) * EunloadP;
			}
		}
		else //sigP < 0
		{
			double zeroSigEps = epsP - sigP / E0n;
			double dE = 0.5 * sigP * (zeroSigEps - epsP);
			EnergyP += dE;
			if (EnergyP < 0) EnergyP = 0.;
			ExcurEnergy += dE;
			if (ExcurEnergy < 0) ExcurEnergy = 0.;
			if (branch == failing)
			{
				FydN = FrN;
				FcN = FrN;
				EshN = (FcN - FydN) / (nd2 - nd1);
			}
			else if (branch == failed)
			{
				FydN = 0;
				FcN = 0;
				EshN = (FcN - FydN) / (nd2 - nd1);
			}
			else {
				double beta = 0;
				if (EnergyP > FailEnergS)
					beta = 1;
				else
					beta = pow(ExcurEnergy / (FailEnergS - EnergyP), cS);
				if (beta > 0.999 || beta < 0)
				{
					opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Complete Strength loss\n" << endln;
					beta = 0.999;
				}
				FydN = (1. - beta) * FydN + beta * nf3 / nf1 * FydN;
				FcN = nf2 - nf1 + FydN;
				EshN = (FcN - FydN) / (nd2 - nd1);
				//FrN = FcN - nf2 + nf3;
				//if (FrN > 0)
				//	FrN = 0;
				if (EnergyP > FailEnergUnloadE)
					beta = 1;
				else
					beta = pow(ExcurEnergy / (FailEnergUnloadE - EnergyP), cUnloadE);
				if (beta > 0.9999)
				{
					opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Complete Unloading Stiffness loss\n" << endln;
					beta = 0.9999;
				}
				EunloadN = (1. - beta) * EunloadN;
			}
		}
		ExcurEnergy = 0.0;
	}
	else
	{
		double dE = 0.5 * (sig + sigP) * (eps - epsP);
		ExcurEnergy += dE;
		EnergyP += dE;
	}
}

int SmoothIMK::setTrialStrain(double trialStrain, double strainRate)
{
	revertToLastCommit();
	eps = trialStrain + (sigini > 0 ? sigini / E0p : sigini / E0n);
	const double deps = eps - epsP;
	if (fabs(deps) < 1.e-20)
		return 0;
	if (!initiated) // handle start in negative direction
	{
		isPosDir = (deps > 0);
		if (gap > 0.001 * pd1_)
		{
			onEnvelope = onEnvelopeP = false;
			branch = gapping;
			if (isPosDir)
			{
				//assuming positive initial loading
				pd1 = pd1_ + gap;
				pd2 = pd2_ + gap;
				pd3 = pd3_ + gap;
				pdu = pdu_ + gap;
				epsmax = pd1;
				e = 0.001 * pf1 / gap;
			}
			else {
				nd1 = nd1_ - gap;
				nd2 = nd2_ - gap;
				nd3 = nd3_ - gap;
				ndu = ndu_ - gap;
				epsmin = nd1;
				e = 0.001 * nf1 / gap;
			}
		}
		initiated = true;
		updateAsymptote();
	}
	if (isPosDir && deps < 0.0)
	{
		isPosDir = false;
		if (sigP > 0)
		{
			epsPl = epsP - sigP / EunloadP;
		}
		if (epsP > epsmax)
			epsmax = epsP;
		onEnvelope = false;
		changeBranch(sigP < 0);
		updateAsymptote();
	}
	else if (deps > 0.0 && !isPosDir) {
		isPosDir = true;
		if (sigP < 0)
		{
			epsPl = epsP - sigP / EunloadN;
		}
		if (epsP < epsmin)
			epsmin = epsP;
		onEnvelope = false;
		changeBranch(sigP > 0);
		updateAsymptote();
	}
	else if ((isPosDir && eps > epsLimit) || (!isPosDir && eps < epsLimit))
	{
		onEnvelope = true;
		changeBranch(false);
		updateAsymptote();
	}
	double epsrat = (eps - epsr) / (epss0 - epsr);
	double dum1 = 1.0 + pow(fabs(epsrat), R0);
	double dum2 = pow(dum1, (1 / R0));

	sig = slopeRat * epsrat + (1.0 - slopeRat) * epsrat / dum2;
	sig = sig * (sigs0 - sigr) + sigr;

	e = slopeRat + (1.0 - slopeRat) / (dum1 * dum2);
	e *= (sigs0 - sigr) / (epss0 - epsr);
	return 0;
}

void SmoothIMK::changeBranch(bool isReturning)
{
	if (onEnvelope)
		branch = nextBranch(branch);
	else if (gap > 0.001 * pd1 && !isReturning)
		branch = pinching;
	else if (isPosDir)
		branch = cyclicRule == 1 ? precap : (cyclicRule == 2 ? peakOriented : pinching);
	//branch = cyclicRule == 1 ? precap : (cyclicRule == 2 ? peakOriented : (epsmax > pd1 ? pinching : precap));
	else
		branch = cyclicRule == 1 ? precap : (cyclicRule == 2 ? peakOriented : pinching);
	//branch = cyclicRule == 1 ? precap : (cyclicRule == 2 ? peakOriented : (epsmin < nd1 ? pinching : precap));
}

void SmoothIMK::updateAsymptote()
{
	switch (cyclicRule)
	{
	case 1:
		bilinAsymptote();
		return;
	case 2:
		peakOrientedAsymptote();
		return;
	case 3:
		pinchedAsymptote();
		return;
		//case 4:
		//	gappedBilinAsymptote();
		//	return;
	default:
		opserr << "SmoothIMK::ERROR Unrecognized cyclic rule: " << cyclicRule << endln;
		exit(-1);

	}
}

void SmoothIMK::bilinAsymptote()
{
	//updates: epsr, sigr, epss0, sigs0, epsLimit, slopeRat, R0, stat
	const double& Esh = isPosDir ? EshP : EshN;
	const double& Fy = isPosDir ? FydP : FydN;
	const double& Fc = isPosDir ? FcP : FcN;
	const double& Fr = isPosDir ? FrP : FrN;
	const double& dy = isPosDir ? pd1 : nd1;
	const double& dc = isPosDir ? pd2 : nd2;
	const double& dr = isPosDir ? pd3 : nd3;
	const double& du = isPosDir ? pdu : ndu;
	const double& E1 = isPosDir ? E0p : E0n;
	const double& E2 = isPosDir ? EunloadN : EunloadP;
	const double& gapD = isPosDir ? gap : -gap;
	//const double& E2 = branchP == peakOriented ? eP : (isPosDir ? E0n : E0p);
	epsr = epsP;
	sigr = sigP;
	double Ec = (Fr - Fc) / (dr - dc);
	double k2 = 0;
	double Dy = 0;
	switch (branch)
	{
	case gapping:
		epss0 = epsPl + gapD;
		sigs0 = 0.001 * Fy;
		Dy = (Fy - Esh * dy - sigs0 + E2 * epss0) / (E2 - Esh);
		epsLimit = (epss0 + Dy) / 2;
		k2 = E1;
		break;
	case precap:
		epss0 = (Fy - Esh * dy - sigr + E2 * epsr) / (E2 - Esh);
		if ((!isPosDir && (epss0 > epsr || epss0 < dc)) || (isPosDir && (epss0 < epsr || epss0 > dc)))
		{
			//we should skip this branch
			branch = postcap;
		}
		else
		{
			sigs0 = Fy + Esh * (epss0 - dy);
			epsLimit = (dy + dc) / 2;
			k2 = Esh;
			break;
		}
	case postcap:
		epss0 = dc;
		sigs0 = Fc;
		epsLimit = (dc + dr) / 2;
		k2 = Ec;
		break;
	case residual:
		epss0 = dr;
		sigs0 = Fr;
		epsLimit = (dr + du) / 2;
		k2 = 0;
		break;
	case failing:
		epss0 = du;
		sigs0 = Fr;
		epsLimit = (3 * du - dr) / 2;
		Dy = du - Fr / Ec;
		if ((isPosDir && Dy < epsLimit) || (!isPosDir && Dy > epsLimit))
			epsLimit = Dy;
		k2 = Ec;
		break;
	case failed:
		epss0 = 1000*du;
		sigs0 = 0;
		epsLimit = 1000 * dy;
		k2 = 0;
		break;
	case pinching:
		sigs0 = 0.001 * Fy;
		epss0 = epsPl;
		epsLimit = epss0 + gapD / 2;
		k2 = 0.001 * E1;
		break;
	}
	double k1 = (sigs0 - sigr) / (epss0 - epsr);
	slopeRat = k2 / k1;
	computeR0(k1, k2, E1, dy);
}

void SmoothIMK::peakOrientedAsymptote()
{
	//updates: epsr, sigr, epss0, sigs0, epsLimit, slopeRat, R0, stat
	if (onEnvelope)
		return bilinAsymptote();
	//const double& E2 = isPosDir ? E0n : E0p;
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	const double& dy = isPosDir ? pd1 : nd1;
	const double& E1 = isPosDir ? E0p : E0n;
	const double& sigPenetFac = isPosDir ? sigPenetFacP : sigPenetFacN;
	const double& Eunload = isPosDir ? EunloadN : EunloadP;
	const double& Esh = isPosDir ? EshP : EshN;
	const double& Fy = isPosDir ? FydP : FydN;
	const double& gapD = isPosDir ? gap : -gap;
	const double& E2 = isPosDir ? E0n : E0p;
	//peakOriented branch
	double tmp;
	double k2;
	epsr = epsP;
	sigr = sigP;
	if (branch == gapping)
	{
		epss0 = epsPl + gapD;
		sigs0 = 0.001 * Fy;
		epsLimit = (epss0 + epsPeak) / 2;
		double sigmax;
		ebranch targBranch; // not used
		double tmp;
		getEnvelope(epsPeak, sigmax, targBranch, tmp, tmp);
		k2 = (sigmax - sigs0) / (epsPeak - epss0);
	}
	else if (branch == pinching) // for gapping
	{
		sigs0 = 0.001 * Fy;
		epss0 = epsPl;
		epsLimit = epss0 + gapD / 2;
		k2 = 0.001 * E1;
	}
	else if (branch == peakOriented)
	{
		double sigmax;
		ebranch targBranch; // not used
		double tmp;
		getEnvelope(epsPeak, sigmax, targBranch, tmp, tmp);
		sigs0 = sigPenetFac * sigmax;
		epss0 = epsPl + sigs0 / Eunload;
		epsLimit = (epss0 + epsPeak) / 2;
		k2 = (sigmax - sigs0) / (epsPeak - epss0);
	}
	else
	{
		epss0 = epsPeak;
		getEnvelope(epss0, sigs0, branch, k2, epsLimit);
	}
	double k1 = (sigs0 - sigr) / (epss0 - epsr);
	slopeRat = k2 / k1;
	computeR0(k1, k2, E1, dy);
}

void SmoothIMK::pinchedAsymptote()
{
	if (branch == gapping || branch == precap || onEnvelope)
		return bilinAsymptote();
	const double& E2 = isPosDir ? E0n : E0p;
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	const double& dy = isPosDir ? pd1 : nd1;
	const double& E1 = isPosDir ? E0p : E0n;
	const double pinchX = isPosDir ? pinchXPos : pinchXNeg;
	const double pinchY = isPosDir ? pinchYPos : pinchYNeg;
	const double& sigPenetFac = isPosDir ? sigPenetFacP : sigPenetFacN;
	const double& Eunload = isPosDir ? EunloadN : EunloadP;
	const double& gapD = isPosDir ? gap : -gap;
	//peakOriented branch
	double tmp;
	double k2;
	epsr = epsP;
	sigr = sigP;
	if (branch == pinching)
	{
		double sigmax;
		ebranch targBranch; // not used
		double tmp;
		getEnvelope(epsPeak, sigmax, targBranch, tmp, tmp);
		sigs0 = sigPenetFac * sigmax;
		epss0 = epsPl + sigs0 / Eunload;
		double x = epsPl + pinchX * (epsPeak - epsPl) + gapD;
		epsLimit = (epss0 + x) / 2;
		double y = pinchY * sigmax;
		k2 = (y - sigs0) / (x - epss0);
	}
	else if (branch == peakOriented)
	{
		epss0 = epsPl + pinchX * (epsPeak - epsPl) + gapD;
		epsLimit = (epss0 + epsPeak) / 2;
		double sigmax;
		ebranch targBranch; // not used
		double tmp;
		getEnvelope(epsPeak, sigmax, targBranch, tmp, tmp);
		sigs0 = pinchY * sigmax;
		k2 = (sigmax - sigs0) / (epsPeak - epss0);
	}
	else
	{
		epss0 = epsPeak;
		getEnvelope(epss0, sigs0, branch, k2, epsLimit);
	}
	double k1 = (sigs0 - sigr) / (epss0 - epsr);
	slopeRat = k2 / k1;
	computeR0(k1, k2, E1, dy);
}

void SmoothIMK::gappedBilinAsymptote()
{
	if (onEnvelope)
		return bilinAsymptote();
	const double& Esh = isPosDir ? EshP : EshN;
	const double& Fy = isPosDir ? FydP : FydN;
	const double& gapD = isPosDir ? gap : -gap;
	const double& E2 = isPosDir ? E0n : E0p;
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	const double& dy = isPosDir ? pd1 : nd1;
	const double& E1 = isPosDir ? E0p : E0n;
	const double pinchX = isPosDir ? pinchXPos : pinchXNeg;
	const double pinchY = isPosDir ? pinchYPos : pinchYNeg;
	const double& sigPenetFac = isPosDir ? sigPenetFacP : sigPenetFacN;
	const double& Eunload = isPosDir ? EunloadN : EunloadP;
	//peakOriented branch
	double tmp;
	double k2;
	epsr = epsP;
	sigr = sigP;
	if (branch == gapping)
	{
		epss0 = epsPl + gapD;
		sigs0 = 0.001 * Fy;
		double Dy = (Fy - Esh * dy - sigs0 + E2 * epss0) / (E2 - Esh);
		epsLimit = (epss0 + Dy) / 2;
		k2 = E1;
	}
	if (branch == pinching)
	{
		sigs0 = 0.001 * Fy;
		epss0 = epsPl;
		epsLimit = epss0 + gapD / 2;
		k2 = 0.001 * E1;
	}
	double k1 = (sigs0 - sigr) / (epss0 - epsr);
	slopeRat = k2 / k1;
	computeR0(k1, k2, E1, dy);
}

SmoothIMK::ebranch SmoothIMK::nextBranch(ebranch branch)
{
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	switch (branch)
	{
	case gapping:
		onEnvelope = true;
		if (cyclicRule == 1)
			return precap;
		else //if (cyclicRule == 2)
			return peakOriented;
	case precap:
		onEnvelope = true;
		return postcap;
	case postcap:
		onEnvelope = true;
		return residual;
	case residual:
		onEnvelope = true;
		if (isPosDir)
		{
			if (FrP == 0)
				return failed;
			else
				return failing;
		}
		else
		{
			if (FrN == 0)
				return failed;
			else
				return failing;
		}
	case failing:
		onEnvelope = true;
		return failed;
	case failed:
		onEnvelope = true;
		return failed;
	case peakOriented:
		onEnvelope = false;
		return peakPassing;
	case pinching:
		onEnvelope = false;
		if (cyclicRule != 3)
			return gapping;
		return peakOriented;
	}
}

void SmoothIMK::getEnvelope(double eps, double& targStress, SmoothIMK::ebranch& targBranch, double& k, double& limitEps)
{
	const double& Esh = eps > 0 ? EshP : EshN;
	const double& Fy = eps > 0 ? FydP : -FydN;
	const double& Fc = eps > 0 ? FcP : -FcN;
	const double& Fr = eps > 0 ? FrP : -FrN;
	const double& dy = eps > 0 ? pd1 : -nd1;
	const double& dc = eps > 0 ? pd2 : -nd2;
	const double& dr = eps > 0 ? pd3 : -nd3;
	const double& du = eps > 0 ? pdu : -ndu;
	const double& E = eps > 0 ? E0p : E0n;
	int sgn = eps > 0 ? 1 : -1;
	eps = fabs(eps);
	if (eps < dy)
	{
		targStress = sgn * E * eps;
		targBranch = precap;
		k = E;
		limitEps = sgn * (dy + dc) / 2.;
		return;
	}
	else if (eps < dc)
	{
		targStress = sgn * (Fy + Esh * (eps - dy));
		k = Esh;
		targBranch = precap;
		limitEps = sgn * (eps + dc) / 2;
		return;
	}
	else if (eps < dr)
	{
		k = (Fr - Fc) / (dr - dc);
		targStress = sgn * (Fc + k * (eps - dc));
		targBranch = postcap;
		limitEps = sgn * (eps + dr) / 2;
		return;
	}
	else if (eps < du)
	{
		targStress = sgn * Fr;
		k = 0;
		targBranch = residual;
		limitEps = sgn * (eps + du) / 2;
		return;
	}
	else // >= du
	{
		targStress = 0;
		targBranch = failed;
		k = 0;
		limitEps = sgn * 1000 * dy;
	}

}

void SmoothIMK::computeR0(double k1, double k2, double E1, double dy)
{
	R0 = r0;
	if (r1 != 0.0)
	{
		double xi_1 = fabs(k1 - k2) / E1;
		double xi_2 = (epss0 - epsr) / dy;
		R0 += r1 * xi_1 + r2 * xi_2;
		R0 = std::min(17.0, R0);
		R0 = std::max(1.0, R0);
		//if (print)
		//	opserr << epsP << " " << xi_1 << " " << xi_2 << " " << R0 << endln;
	}
}
