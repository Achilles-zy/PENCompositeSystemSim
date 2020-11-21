#include "PENDetectorConstruction.hh"

#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Hype.hh"
#include "G4EllipticalCone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4String.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4NistManager.hh"
#include "G4Cerenkov.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4MultiUnion.hh"

#include "PENMaterials.hh"
#include <G4VisAttributes.hh>
#include <iostream>
#include <fstream>
#include <iterator>

#include "CADMesh.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

//#include "TMath.h"

using namespace std;
PENDetectorConstruction::PENDetectorConstruction():
	G4VUserDetectorConstruction(),
	physEnv(nullptr),
	physBulk(nullptr),
	physPENShell(nullptr),
    physSiPM0(nullptr),
    physSiPM1(nullptr),
    physSiPM2(nullptr),
    physSiPM3(nullptr),
    physSiPM4(nullptr),
    physSiPM5(nullptr),
    physSiPM6(nullptr),
    physSiPM7(nullptr),
    physSiPM8(nullptr),
    physSiPM9(nullptr),
	physSiPMArray0(nullptr),
	physSiPMArray1(nullptr),
	physSiPMArray2(nullptr),
	physSiPMArray3(nullptr),
	physOuterReflector(nullptr),
	physInnerReflector(nullptr),
	physWire(nullptr),
	physSArBrick(nullptr),
	//logicPENShell(nullptr),
	solidSideSiPM(nullptr)
{
	fDetectorMessenger = new PENDetectorMessenger(this);
	matconstructor = new PENMaterials;
	MPT_PEN = new G4MaterialPropertiesTable();
	AbsorptionLength = 1.5;//value at 400 nm
	fRES = 1.0;
	fLY = 3500. / MeV;
	fABSFile = "PEN_ABS";
	fConfine = "PENShell";
	fWireType = "A1";
	fReflectorType = "ESR";
	fMode = "SArUnit";
	fWirePos = G4ThreeVector();
	fWireRadius = 0.7 * mm;
	fWireLength = 20 * cm;
	fWireCentDist = 5 * cm;
	pmtReflectivity = 0.50;
	fPENShellLength = 10 * cm;
	fPENShellRadius = 20 * cm;
	fPENPropertiesID = 1;
	fBEGeRadius = 40 * mm;
	fBEGeHeight = 40 * mm;
	fSArBrickHeight = 150 * cm;
	fSArBrickRadius = 40 * cm;
	absFactor = 1.5;
	fDetMat = matEnGe;
	G4cout << "Start Construction" << G4endl;
	DefineMat();
	fTargetMaterial = G4Material::GetMaterial("PVT_structure");
	fGlassMaterialPMT = G4Material::GetMaterial("BorosilicateGlass");
	ifOuterReflector = false;
	ifInnerReflector = false;
	CheckOverlaps = true;
}

PENDetectorConstruction::~PENDetectorConstruction()
{
}

void PENDetectorConstruction::SetWireType(G4String type) {
	fWireType = type;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void PENDetectorConstruction::SetConfine(G4String confine) {
	fConfine = confine;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void PENDetectorConstruction::SetRunInfo(G4String runinfo) {
	fRunInfo = runinfo;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void PENDetectorConstruction::SetMode(G4String mode) {
	fMode = mode;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void PENDetectorConstruction::SetPENPropertiesID(G4int nb) {
	fPENPropertiesID = nb;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void PENDetectorConstruction::SetOuterReflector(G4bool ref) {
	ifOuterReflector = ref;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void PENDetectorConstruction::SetInnerReflector(G4bool ref) {
	ifInnerReflector = ref;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void PENDetectorConstruction::SetReflectorType(G4String type) {
	fReflectorType = type;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void PENDetectorConstruction::DefineMat() 
{
	matconstructor->Construct();
	// ============================================================= Materials =============================================================
  //materialConstruction = new PenMaterials;
	matAir = G4Material::GetMaterial("Air");
	matBialkali = G4Material::GetMaterial("Bialkali");
	fGlass = G4Material::GetMaterial("BorosilicateGlass");
	fPOM = G4Material::GetMaterial("POM");
	fABS = G4Material::GetMaterial("ABS");
	matPEN = G4Material::GetMaterial("PEN");
	matSi = G4Material::GetMaterial("G4_Si");
	matCu = G4Material::GetMaterial("G4_Cu");
	matTriggerFoilEJ212 = G4Material::GetMaterial("EJ212");
	Pstyrene = G4Material::GetMaterial("Polystyrene");
	matPMMA = G4Material::GetMaterial("PMMA");
	fVacuum = G4Material::GetMaterial("Vacuum");
	matGreaseEJ550 = G4Material::GetMaterial("Grease");
	matTeflon = G4Material::GetMaterial("G4_TEFLON");
	matVikuiti = G4Material::GetMaterial("Vikuiti");
	matTitanium = G4Material::GetMaterial("titanium");
	matPolyethylene = G4Material::GetMaterial("G4_POLYETHYLENE");
	matEnGe = G4Material::GetMaterial("EnGe");
	matNaGe = G4Material::GetMaterial("NaGe");
	matLN2 = G4Material::GetMaterial("G4_lN2");
	matLAr = G4Material::GetMaterial("G4_lAr");
	matGAGG = G4Material::GetMaterial("GAGG");
	matPTFE = G4Material::GetMaterial("PTFE");

	G4cout << " materials ok " << G4endl;

	G4double wavelength;
	char filler;
	G4double varAbsorLength;
	G4double emission;
	G4double rindex;

	G4double wlPhotonEnergy[102] = { 0 };
	G4double ABSORPTION_PEN[102] = { 0 };
	G4double RINDEX_PEN[102] = { 0 };

	G4int absEntries = 0;

	ifstream ReadAbs;

	G4String absFile = "../input_files/" + fABSFile + ".csv";
	ReadAbs.open(absFile);
	if (ReadAbs.is_open())
	{
		while (!ReadAbs.eof())
		{
			ReadAbs >> wavelength >> filler >> varAbsorLength >> filler >> emission >> filler >> rindex;
			if (ReadAbs.eof()) {
				break;
			}
			wlPhotonEnergy[absEntries] = (1240. / wavelength) * eV;
			ABSORPTION_PEN[absEntries] = (varAbsorLength) * mm;
			RINDEX_PEN[absEntries] = rindex;
			absEntries++;
		}
	}

	else G4cout << "Error opening file: " << absFile << G4endl;
	ReadAbs.close();
	absEntries--;

	const G4int nEntries1 = sizeof(wlPhotonEnergy) / sizeof(G4double);
	assert(sizeof(RINDEX_PEN) == sizeof(wlPhotonEnergy));
	assert(sizeof(ABSORPTION_PEN) == sizeof(wlPhotonEnergy));
	//assert(sizeof(EMISSION_PEN) == sizeof(wlPhotonEnergy));

	MPT_PEN = new G4MaterialPropertiesTable();

	MPT_PEN->AddProperty("RINDEX", wlPhotonEnergy, RINDEX_PEN, nEntries1)->SetSpline(true);
	MPT_PEN->AddProperty("ABSLENGTH", wlPhotonEnergy, ABSORPTION_PEN, nEntries1)->SetSpline(true); // *

	// Read primary emission spectrum from PEN
	// Measurements from MPP Munich
	G4double pWavelength;
	G4String  Scint_file = "../properties/PEN_EM_SPECTRUM.dat";
	std::ifstream ReadScint2(Scint_file), ReadScintPEN;
	//count number of entries
	ReadScint2.unsetf(std::ios_base::skipws);
	//unsigned line_count = std::count(
	int line_count = std::count(
		std::istream_iterator<char>(ReadScint2),
		std::istream_iterator<char>(),
		'\n');
	std::cout << "Lines: " << line_count << "\n";
	ReadScint2.close();
	G4double PEN_EMISSION[500];
	G4double PEN_WL_ENERGY[500];
	G4int nEntriesPEN = 0;
	ReadScintPEN.open(Scint_file);
	if (ReadScintPEN.is_open()) {
		while (!ReadScintPEN.eof()) {

			ReadScintPEN >> pWavelength >> PEN_EMISSION[nEntriesPEN];
			if (ReadScintPEN.eof()) {
				break;
			}
			PEN_WL_ENERGY[nEntriesPEN] = (1240. / pWavelength) * eV;//convert wavelength to eV
		//G4cout<<nEntriesPEN<<" wl "<<PEN_WL_ENERGY[nEntriesPEN]<<" "<<PEN_EMISSION[nEntriesPEN]<<G4endl;
			nEntriesPEN++;
			if (nEntriesPEN > (line_count - 1)) { G4cout << " entries completed " << G4endl; break; }
		}
	}
	else
		G4cout << "Error opening file: " << Scint_file << G4endl;
	ReadScintPEN.close();
	G4cout << " nEntriesPEN " << nEntriesPEN << G4endl;

	MPT_PEN->AddProperty("FASTCOMPONENT", PEN_WL_ENERGY, PEN_EMISSION, line_count)->SetSpline(true);
	MPT_PEN->AddProperty("SLOWCOMPONENT", PEN_WL_ENERGY, PEN_EMISSION, line_count)->SetSpline(true);

	MPT_PEN->AddConstProperty("SCINTILLATIONYIELD", fLY); // * 2.5 * PEN = PS, 10*PEN=PS
	MPT_PEN->AddConstProperty("RESOLUTIONSCALE", fRES); // * 1, 4, 8
	MPT_PEN->AddConstProperty("FASTTIMECONSTANT", 5.198 * ns);
	MPT_PEN->AddConstProperty("SLOWTIMECONSTANT", 24.336 * ns);
	MPT_PEN->AddConstProperty("YIELDRATIO", 1.);

	G4cout << "PEN Properties:" << G4endl;
	G4cout << "AbsFactor =" << absFactor << G4endl;
	G4cout << "LY =" << fLY << G4endl;
	matPEN->SetMaterialPropertiesTable(MPT_PEN);
	//pvt_structure->SetMaterialPropertiesTable(MPT_PEN);


	G4cout << " pen ok " << G4endl;


	G4double rindexEnergy[500] = { 0 };
	G4double scintIndex[500] = { 0 };

	G4int rindexEntries = 0;
	ifstream ReadRindex;

	G4String rindex_file = "../input_files/rindexScint.txt";
	ReadRindex.open(rindex_file);

	if (ReadRindex.is_open())
	{
		while (!ReadRindex.eof())
		{

			ReadRindex >> wavelength >> filler >> scintIndex[rindexEntries];
			if (ReadRindex.eof()) {
				break;
			}
			rindexEnergy[rindexEntries] = (1240. / wavelength) * eV;
			rindexEntries++;
		}
	}
	else G4cout << "Error opening file: " << rindex_file << G4endl;
	ReadRindex.close();
	rindexEntries--;

	G4double scintEnergy[501] = { 0 };
	G4double scintEmit[501] = { 0 };
	G4double scintEmitSlow[501] = { 0 };

	G4int scintEntries = 0;
	ifstream ReadScint;

	Scint_file = "../input_files/pTP_emission.txt";
	ReadScint.open(Scint_file);

	if (ReadScint.is_open())
	{
		while (!ReadScint.eof())
		{

			ReadScint >> wavelength >> filler >> scintEmit[scintEntries];
			if (ReadScint.eof()) {
				break;
			}
			//convert wavelength to eV:
			scintEnergy[scintEntries] = (1240. / wavelength) * eV;
			scintEmitSlow[scintEntries] = scintEmit[scintEntries];
			scintEntries++;
		}
	}
	else G4cout << "Error opening file: " << Scint_file << G4endl;
	ReadScint.close();
	scintEntries--;

	G4int absorbEntries = 0;
	G4double varAbsorbLength;
	G4double absorbEnergy[501] = { 0 };
	G4double Absorb[501] = { 0 };

	ifstream ReadAbsorb;
	G4String ReadAbsorbLength = "../input_files/PlasticBulkAbsorb2.cfg";

	ReadAbsorb.open(ReadAbsorbLength);
	if (ReadAbsorb.is_open())
	{
		while (!ReadAbsorb.eof())
		{

			ReadAbsorb >> wavelength >> filler >> varAbsorbLength;
			if (ReadAbsorb.eof()) {
				break;
			}
			absorbEnergy[absorbEntries] = (1240 / wavelength) * eV;
			Absorb[absorbEntries] = (varAbsorbLength)*m;
			absorbEntries++;
		}
	}
	else G4cout << "Error opening file: " << ReadAbsorbLength << G4endl;
	ReadAbsorb.close();
	absorbEntries--;

	G4double wlsEnergy[501] = { 0 };
	G4double wlsEmit[501] = { 0 };

	G4int wlsScintEntries = 0;
	ifstream ReadWLSScint;

	G4String wls_Scint_file = "../input_files/full_popop_emission.cfg";
	ReadWLSScint.open(wls_Scint_file);

	if (ReadWLSScint.is_open())
	{
		while (!ReadWLSScint.eof())
		{

			ReadWLSScint >> wavelength >> filler >> wlsEmit[500 - wlsScintEntries];
			if (ReadWLSScint.eof()) {
				break;
			}
			//convert wavelength to eV:
			wlsEnergy[500 - wlsScintEntries] = (1240 / wavelength) * eV;
			wlsScintEntries++;
		}
	}
	else G4cout << "Error opening file: " << wls_Scint_file << G4endl;
	ReadWLSScint.close();
	wlsScintEntries--;

	G4int wlsAbsorbEntries = 0;
	G4double wlsAbsorbEnergy[501] = { 0 };
	G4double wlsAbsorb[501] = { 0 };

	ifstream ReadWLSAbsorb;
	G4String ReadWLSAbsorbLength = "../input_files/scintAbsLen.txt";

	ReadWLSAbsorb.open(ReadWLSAbsorbLength);
	if (ReadWLSAbsorb.is_open())
	{
		while (!ReadWLSAbsorb.eof())
		{
			ReadWLSAbsorb >> wavelength >> filler >> varAbsorbLength;
			if (ReadWLSAbsorb.eof()) {
				break;
			}
			wlsAbsorbEnergy[wlsAbsorbEntries] = (1240. / wavelength) * eV;
			wlsAbsorb[wlsAbsorbEntries] = varAbsorbLength * mm;
			wlsAbsorbEntries++;
		}
	}
	else G4cout << "Error opening file: " << ReadWLSAbsorbLength << G4endl;
	ReadWLSAbsorb.close();
	wlsAbsorbEntries--;

	G4MaterialPropertiesTable* MPT_FoilEJ212 = new G4MaterialPropertiesTable();

	MPT_FoilEJ212->AddProperty("WLSABSLENGTH", wlsAbsorbEnergy, wlsAbsorb, wlsAbsorbEntries);
	MPT_FoilEJ212->AddProperty("WLSCOMPONENT", wlsEnergy, wlsEmit, wlsScintEntries);
	MPT_FoilEJ212->AddConstProperty("WLSTIMECONSTANT", 12 * ns);

	MPT_FoilEJ212->AddProperty("RINDEX", rindexEnergy, scintIndex, rindexEntries);
	MPT_FoilEJ212->AddProperty("ABSLENGTH", absorbEnergy, Absorb, absorbEntries);
	MPT_FoilEJ212->AddProperty("FASTCOMPONENT", scintEnergy, scintEmit, scintEntries);
	MPT_FoilEJ212->AddProperty("SLOWCOMPONENT", scintEnergy, scintEmitSlow, scintEntries);

	//MPT_FoilEJ212->AddConstProperty("SCINTILLATIONYIELD",11520./MeV);
	MPT_FoilEJ212->AddConstProperty("SCINTILLATIONYIELD", 10. / MeV);//set low LY to make it faster, intead use Edep for coincidences
	MPT_FoilEJ212->AddConstProperty("RESOLUTIONSCALE", 4.0);
	MPT_FoilEJ212->AddConstProperty("FASTTIMECONSTANT", 2.1 * ns);
	MPT_FoilEJ212->AddConstProperty("SLOWTIMECONSTANT", 14.2 * ns);
	MPT_FoilEJ212->AddConstProperty("YIELDRATIO", 1.0);

	matTriggerFoilEJ212->SetMaterialPropertiesTable(MPT_FoilEJ212);

	G4cout << " EJ212 ok " << G4endl;

	G4double refractive_index[] = { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49 };
	G4double absPMMA[] = { 5 * m, 5 * m, 5 * m, 5 * m, 5 * m, 5 * m };
	G4double reflPMMA[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
	G4double energyPMMA[] = { 2.18 * eV, 2.48 * eV, 2.58 * eV, 2.68 * eV, 2.78 * eV, 4.1 * eV };
	const G4int nEntries3 = sizeof(energyPMMA) / sizeof(G4double);

	G4MaterialPropertiesTable* MPT_PMMA = new G4MaterialPropertiesTable();
	MPT_PMMA->AddProperty("RINDEX", energyPMMA, refractive_index, nEntries3);
	MPT_PMMA->AddProperty("ABSLENGTH", energyPMMA, absPMMA, nEntries3)->SetSpline(true);
	MPT_PMMA->AddProperty("REFLECTIVITY", energyPMMA, reflPMMA, nEntries3)->SetSpline(true);
	matPMMA->SetMaterialPropertiesTable(MPT_PMMA);

}

void PENDetectorConstruction::SetABS(G4double value) {
	AbsorptionLength = value;
	//read file and add the value given by the user
	G4double wavelength;
	char filler;
	G4double varAbsorLength;
	G4double emission;
	G4double rindex;

	G4double wlPhotonEnergy[102] = { 0 };
	G4double ABSORPTION_PEN[102] = { 0 };

	G4int absEntries = 0;
	ifstream ReadAbs;

	G4String absFile = "../input_files/" + fABSFile + ".csv";
	ReadAbs.open(absFile);
	if (ReadAbs.is_open())
	{
		while (!ReadAbs.eof())
		{
			ReadAbs >> wavelength >> filler >> varAbsorLength >> filler >> emission >> filler >> rindex;
			if (ReadAbs.eof()) {
				break;
			}
			wlPhotonEnergy[absEntries] = (1240. / wavelength) * eV;
			ABSORPTION_PEN[absEntries] = (varAbsorLength * AbsorptionLength) * mm; //use measured value of attenuation to constrain curve and then change values multiplying the curve for a given factor
			absEntries++;
		}
	}

	else G4cout << "Error opening file: " << absFile << G4endl;
	ReadAbs.close();
	absEntries--;

	const G4int nEntries1 = sizeof(wlPhotonEnergy) / sizeof(G4double);
	assert(sizeof(ABSORPTION_PEN) == sizeof(wlPhotonEnergy));
	MPT_PEN->AddProperty("ABSLENGTH", wlPhotonEnergy, ABSORPTION_PEN, nEntries1)->SetSpline(true); // *
	//G4RunManager::GetRunManager()->PhysicsHasBeenModified();
#ifdef G4MULTITHREADED
	G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
#else
	G4RunManager::GetRunManager()->PhysicsHasBeenModified();
#endif
}

void PENDetectorConstruction::SetLY(G4double ly) {

	MPT_PEN->AddConstProperty("SCINTILLATIONYIELD", ly); // * 2.5 * PEN = PS, 10*PEN=PS
//G4RunManager::GetRunManager()->PhysicsHasBeenModified();
#ifdef G4MULTITHREADED
	G4MTRunManager::GetRunManager()->PhysicsHasBeenModified();
#else
	G4RunManager::GetRunManager()->PhysicsHasBeenModified();
#endif
}

////////
//
////////


G4LogicalVolume* PENDetectorConstruction::ConstructBEGe() {

	//flatBEGe
	G4bool checkOverlaps = true;
	G4double GeChamfer = 2. * mm;
	//double outerGeRadius = 31.35 * mm;
	G4double outerGeRadius = fBEGeRadius;
	G4double innerGeRadius = 1.50 * mm;
	//double GeHeight1 = 60.80 * mm;
	G4double GeHeight1 = fBEGeHeight - GeChamfer;
	G4double SrcThickness = 0.01 * mm;
	G4double lSmallValue = 0.01 * mm;
	G4double groovedepth = 1.5 * mm;
	G4double grooveradius = 9 * mm;
	G4double groovethickness = 0.001 * mm;
	G4double outerplayerthickness = 1 * um;
	G4double orbradius = 1.89 * mm;
	G4double deadlayerthickness = 1 * mm;

	// G4Height2 is the depth of the hole for pin contact
	double GeHeight2 = 4.0 * mm;
	double GeHeight3 = GeHeight2 + deadlayerthickness;

	G4ThreeVector zTransGe0(0., 0., -GeHeight1 / 2);
	G4ThreeVector zTransGe1(0., 0., -GeHeight1 / 2 + deadlayerthickness);
	G4ThreeVector zTransGe2(0., 0., GeHeight1 / 2 + lSmallValue);
	G4ThreeVector zTransGroove(0., 0., -GeHeight1 / 2);

	auto GeP1 = new G4Tubs("GeP1", 0., outerGeRadius, GeHeight1 / 2, 0., twopi);
	auto GeP2 = new G4Tubs("GeP2", 0., outerGeRadius - GeChamfer, GeChamfer, 0., twopi);
	auto GeP3 = new G4Torus("GeP3", 0., GeChamfer, outerGeRadius - GeChamfer, 0., twopi);
	auto GeGroove = new G4Torus("solidgroove", 0., groovedepth, grooveradius, 0., twopi);

	// total germanium crystal
	auto GeTemp1 = new G4UnionSolid("GeTemp1", GeP2, GeP3);
	auto GeTemp2 = new G4UnionSolid("GeTemp2", GeP1, GeTemp1, 0, zTransGe2);
	//?
	//auto solidtempTotalCrystal = new G4SubtractionSolid("totaltempCrystal", GeTemp2, GeGroove, 0, zTransGroove);
	auto solidTotalCrystal = new G4SubtractionSolid("totalCrystal", GeTemp2, GeGroove, 0, zTransGroove);
	auto logicTotalCrystal = new G4LogicalVolume(solidTotalCrystal, matEnGe, "logicTotalCrystal");

	// bulk
	auto GeP1In = new G4Tubs("GeP1In", 0., outerGeRadius - deadlayerthickness, (GeHeight1 - deadlayerthickness * 2) / 2, 0., twopi);
	auto GeP2In = new G4Tubs("GeP2In", 0., outerGeRadius - deadlayerthickness - GeChamfer, GeChamfer, 0., twopi);
	auto GeP3In = new G4Torus("GeP3In", 0., GeChamfer, outerGeRadius - deadlayerthickness - GeChamfer, 0., twopi);
	auto GeP4In = new G4Tubs("GeP4In", 0., grooveradius, deadlayerthickness / 2, 0., twopi);

	auto GeLargerGroove = new G4Torus("solidlargergroove", 0., groovedepth + groovethickness, grooveradius, 0., twopi);
	auto GeOuterpLayer = new G4Tubs("solidouterplayer", 0., grooveradius - groovedepth, outerplayerthickness / 2, 0., twopi);

	auto GeInTemp1 = new G4UnionSolid("GeInTemp1", GeP2In, GeP3In);
	auto GeInTemp2 = new G4UnionSolid("GeInTemp2", GeP1In, GeInTemp1, 0, G4ThreeVector(0., 0., (GeHeight1 - deadlayerthickness * 2) / 2.0 + lSmallValue));
	auto GeInTemp3 = new G4UnionSolid("GeInTemp3", GeInTemp2, GeP4In, 0, G4ThreeVector(0., 0., -(GeHeight1 - deadlayerthickness) / 2));
	G4ThreeVector zbulkTrans(0., 0., -(GeHeight1 - deadlayerthickness * 2 - (GeHeight2 - deadlayerthickness)) / 2);
	G4ThreeVector zbulkTransOuterp(0., 0., -(GeHeight1 - outerplayerthickness) / 2);

	auto GeInTemp4 = new G4SubtractionSolid("GeInTemp4", GeInTemp3, GeLargerGroove, 0, zTransGroove);
	auto solidBulk = new G4SubtractionSolid("Bulk", GeInTemp4, GeOuterpLayer, 0, zbulkTransOuterp);
	auto logicBulk = new G4LogicalVolume(solidBulk, matEnGe, "Bulk");

	//deadlayer
	auto tempdeadlayer = new G4SubtractionSolid("tempdeadlayer", solidTotalCrystal, GeInTemp3, 0, G4ThreeVector(0., 0., 0.));
	auto solidOuterDeadlayer = new G4SubtractionSolid("OuterDeadlayer", tempdeadlayer, GeLargerGroove, 0, zTransGroove);
	auto logicOuterDeadlayer = new G4LogicalVolume(solidOuterDeadlayer, matEnGe, "OuterDeadlayer");

	//groove Layer
	auto GrooveTorus = new G4Torus("solidpLayer", groovedepth, groovedepth + groovethickness, grooveradius, 0., twopi);
	auto GrooveCut1 = new G4Tubs("groovecut1", 0, grooveradius + groovedepth + groovethickness, groovedepth + groovethickness, 0, twopi);
	auto GrooveCut2 = new G4Tubs("groovecut2", 0, grooveradius - groovedepth, outerplayerthickness, 0, twopi);
	auto tempGroove1 = new G4SubtractionSolid("tempgroove1", GrooveTorus, GrooveCut1, 0, G4ThreeVector(0., 0., -(groovedepth + groovethickness)));
	auto GrooveLayer = new G4SubtractionSolid("GrooveLayer", tempGroove1, GrooveCut2, 0, G4ThreeVector(0., 0., 0.));
	auto logicGrooveLayer = new G4LogicalVolume(GrooveLayer, matEnGe, "GrooveLayer");

	//outer pLayer
	auto OuterpLayer = new G4Tubs("solidOuterpLayer", 0, grooveradius - groovedepth, outerplayerthickness / 2, 0, twopi);
	auto logicOuterpLayer = new G4LogicalVolume(OuterpLayer, matEnGe, "OuterpLayer");
	physBulk = new G4PVPlacement(0, G4ThreeVector(), logicBulk, "Bulk", logicTotalCrystal, false, 0, checkOverlaps);
	//new G4PVPlacement(0, G4ThreeVector(), logicOuterDeadlayer, "OuterDeadlayer", logicTotalCrystal, false, 0, checkOverlaps);
	//new G4PVPlacement(0, zbulkTransOuterp, logicOuterpLayer, "OuterpLayer", logicTotalCrystal, false, 0, checkOverlaps);
	//new G4PVPlacement(0, zTransGroove, logicGrooveLayer, "GrooveLayer", logicTotalCrystal, false, 0, checkOverlaps);
	return logicTotalCrystal;
}

G4LogicalVolume* PENDetectorConstruction::ConstructA1(G4double WireLength) {
	//=======================================
	// A1 rho = 2.8457 g/m, 2.80 g/m from ref.
	//=======================================

	G4double ConductorRadius = 0.0762 * mm / 2;
	G4double JacketRadius = 0.254 * mm / 2;
	G4double BraidThickness = 0.02 * mm;
	G4double BraidRadius = JacketRadius + BraidThickness;
	G4double SubWireRadius = 0.5 * mm / 2;
	G4double WireRadius = 1.25 * mm / 2;
	fWireRadius = WireRadius;
	G4double WireJacketThickness = WireRadius - (1 + sqrt(2)) * SubWireRadius;
	G4Material* JacketMat = matPTFE;

	G4Tubs* solidWire = new G4Tubs("solidWire", 0 * mm, WireRadius, WireLength / 2, 0, twopi);
	G4LogicalVolume* logicWire = new G4LogicalVolume(solidWire, JacketMat, "logicWire");

	G4Tubs* solidWireJacket = new G4Tubs("solidWireJacket", WireRadius - WireJacketThickness, WireRadius, WireLength / 2, 0, twopi);
	G4LogicalVolume* logicWireJacket = new G4LogicalVolume(solidWireJacket, JacketMat, "logicWireJacket");
	G4PVPlacement* physWireJacket = new G4PVPlacement(0, G4ThreeVector(), logicWireJacket, "WireJacket", logicWire, false, 0, CheckOverlaps);

	G4Tubs* solidSubWire = new G4Tubs("solidWire", 0 * mm, SubWireRadius, WireLength / 2, 0, twopi);
	G4LogicalVolume* logicSubWire = new G4LogicalVolume(solidSubWire, JacketMat, "logicWire");
	G4PVPlacement* physSubWire_0 = new G4PVPlacement(0, G4ThreeVector(sqrt(2) * SubWireRadius, 0, 0), logicSubWire, "SubWire_0", logicWire, false, 0, CheckOverlaps);
	G4PVPlacement* physSubWire_1 = new G4PVPlacement(0, G4ThreeVector(0, sqrt(2) * SubWireRadius, 0), logicSubWire, "SubWire_1", logicWire, false, 1, CheckOverlaps);
	G4PVPlacement* physSubWire_2 = new G4PVPlacement(0, G4ThreeVector(-sqrt(2) * SubWireRadius, 0, 0), logicSubWire, "SubWire_2", logicWire, false, 2, CheckOverlaps);
	G4PVPlacement* physSubWire_3 = new G4PVPlacement(0, G4ThreeVector(0, -sqrt(2) * SubWireRadius, 0), logicSubWire, "SubWire_3", logicWire, false, 3, CheckOverlaps);

	G4Tubs* solidConductor = new G4Tubs("solidConductor", 0 * mm, ConductorRadius, WireLength / 2, 0, twopi);
	G4LogicalVolume* logicConductor = new G4LogicalVolume(solidConductor, matCu, "logicConductor");
	G4PVPlacement* physConductor = new G4PVPlacement(0, G4ThreeVector(), logicConductor, "Conductor", logicSubWire, false, 0, CheckOverlaps);

	G4Tubs* solidJacket = new G4Tubs("solidJacket", ConductorRadius, JacketRadius, WireLength / 2, 0, twopi);
	G4LogicalVolume* logicJacket = new G4LogicalVolume(solidJacket, JacketMat, "logicJacket");
	G4PVPlacement* physJacket = new G4PVPlacement(0, G4ThreeVector(), logicJacket, "Jacket", logicSubWire, false, 0, CheckOverlaps);

	G4Tubs* solidBraid = new G4Tubs("solidBraid", JacketRadius, BraidRadius, WireLength / 2, 0, twopi);
	G4LogicalVolume* logicBraid = new G4LogicalVolume(solidBraid, matCu, "logicBraid");
	G4PVPlacement* physBraid = new G4PVPlacement(0, G4ThreeVector(), logicBraid, "Braid", logicSubWire, false, 0, CheckOverlaps);

	G4Tubs* solidOuterJacket = new G4Tubs("solidOuterJacket", BraidRadius, SubWireRadius, WireLength / 2, 0, twopi);
	G4LogicalVolume* logicOuterJacket = new G4LogicalVolume(solidOuterJacket, JacketMat, "logicOuterJacket");
	G4PVPlacement* physOuterJacket = new G4PVPlacement(0, G4ThreeVector(), logicOuterJacket, "OuterJacket", logicSubWire, false, 0, CheckOverlaps);
	return logicWire;
}

G4LogicalVolume* PENDetectorConstruction::ConstructA2(G4double WireLength) {
		  //=======================================
	  // A2 rho = 2.0531 g/m, 2 g/m from ref.
	  //=======================================

	  G4double ConductorRadius = 0.0762 * mm / 2;
	  G4double BraidThickness = 0.04 * mm;
	  G4double JacketRadius = 0.68 * mm / 2;
	  G4double BraidRadius = JacketRadius + BraidThickness;
	  G4double WireRadius = 1 * mm / 2;
	  fWireRadius = WireRadius;
	  G4Material* JacketMat = matPTFE;

	  G4Tubs* solidWire = new G4Tubs("solidJacket", 0 * mm, WireRadius, WireLength / 2, 0, twopi);
	  G4LogicalVolume* logicWire = new G4LogicalVolume(solidWire, JacketMat, "logicWire");

	  G4Tubs* solidConductor = new G4Tubs("solidConductor", 0 * mm, ConductorRadius, WireLength / 2, 0, twopi);
	  G4LogicalVolume* logicConductor = new G4LogicalVolume(solidConductor, matCu, "logicConductor");

	  G4Tubs* solidJacket = new G4Tubs("solidJacket", ConductorRadius, JacketRadius, WireLength / 2, 0, twopi);
	  G4LogicalVolume* logicJacket = new G4LogicalVolume(solidJacket, JacketMat, "logicJacket");

	  G4Tubs* solidBraid = new G4Tubs("solidBraid", JacketRadius, BraidRadius, WireLength / 2, 0, twopi);
	  G4LogicalVolume* logicBraid = new G4LogicalVolume(solidBraid, matCu, "logicBraid");

	  G4Tubs* solidOuterJacket = new G4Tubs("solidOuterJacket", BraidRadius, WireRadius, WireLength / 2, 0, twopi);
	  G4LogicalVolume* logicOuterJacket = new G4LogicalVolume(solidOuterJacket, JacketMat, "logicOuterJacket");


	  G4PVPlacement* physConductor = new G4PVPlacement(0, G4ThreeVector(), logicConductor, "Conductor", logicWire, false, 0, CheckOverlaps);
	  G4PVPlacement* physJacket = new G4PVPlacement(0, G4ThreeVector(), logicJacket, "Jacket", logicWire, false, 0, CheckOverlaps);
	  G4PVPlacement* physBraid = new G4PVPlacement(0, G4ThreeVector(), logicBraid, "Braid", logicWire, false, 0, CheckOverlaps);
	  G4PVPlacement* physOuterJacket = new G4PVPlacement(0, G4ThreeVector(), logicOuterJacket, "OuterJacket", logicWire, false, 0, CheckOverlaps);
  
	  return logicWire;
}

G4LogicalVolume* PENDetectorConstruction::ConstructPENShell() {
	auto innermesh = CADMesh::TessellatedMesh::FromSTL("../models/InnerShell-I.stl");
	auto outermesh = CADMesh::TessellatedMesh::FromSTL("../models/OuterShell-II.stl");
	G4VSolid* solidInnerShell = innermesh->GetSolid();
	G4VSolid* solidOuterShell = outermesh->GetSolid();
	//G4ThreeVector position1 = G4ThreeVector(0., -4 * mm, 0.);
	//auto solidPENShell = new G4UnionSolid("solidPENShell", solidInnerShell, solidOuterShell, 0, position1);
	G4ThreeVector position1 = G4ThreeVector(0., 4 * mm, 0.);
	auto solidPENShell = new G4UnionSolid("solidPENShell", solidOuterShell, solidInnerShell, 0, position1);

	/*
	G4MultiUnion* solidPENShell = new G4MultiUnion("solidPENShell");
	// Add the shapes to the structure
	//
	G4RotationMatrix rotm = G4RotationMatrix();
	G4ThreeVector u = G4ThreeVector(1, 0, 0);
	G4ThreeVector v = G4ThreeVector(0, 0, -1);
	G4ThreeVector w = G4ThreeVector(0, 1, 0);
	G4RotationMatrix rotm1 = G4RotationMatrix(u, v, w);
	G4ThreeVector position0 = G4ThreeVector(0., 0., 0.);
	G4ThreeVector position1 = G4ThreeVector(0., 0., 4.1 * mm);
	G4Transform3D tr0 = G4Transform3D(rotm1, position0);
	G4Transform3D tr1 = G4Transform3D(rotm1, position1);

	solidPENShell->AddNode(*solidInnerShell, tr0);
	solidPENShell->AddNode(*solidOuterShell, tr1);
	solidPENShell->Voxelize();
	*/
	auto logicPENShell = new G4LogicalVolume(solidPENShell, matPEN, "logicPENShell");
	auto logicOuterShell = new G4LogicalVolume(solidOuterShell, matPEN, "logicOuterShell");
	auto logicInnerShell = new G4LogicalVolume(solidInnerShell, matPEN, "logicInnerShell");
	return logicPENShell;
}

G4LogicalVolume* PENDetectorConstruction::ConstructOuterReflector() {
	auto outermesh = CADMesh::TessellatedMesh::FromSTL("../models/OuterShell-II.stl");
	G4VSolid* solidOuterShell = outermesh->GetSolid();
	//G4ThreeVector position1 = G4ThreeVector(0., -4 * mm, 0.);
	//auto solidPENShell = new G4UnionSolid("solidPENShell", solidInnerShell, solidOuterShell, 0, position1);

	auto reflectormesh = CADMesh::TessellatedMesh::FromSTL("../models/OuterShell-II-OuterReflector.stl");
	G4VSolid* solidTemp = reflectormesh->GetSolid();
	auto solidReflector = new G4SubtractionSolid("solidReflector", solidTemp, solidOuterShell, 0, G4ThreeVector());
	auto logicReflector = new G4LogicalVolume(solidReflector, matLN2, "logicReflector");
	return logicReflector;

}

G4LogicalVolume* PENDetectorConstruction::ConstructInnerReflector() {
	auto outermesh = CADMesh::TessellatedMesh::FromSTL("../models/OuterShell-II.stl");
	G4VSolid* solidOuterShell = outermesh->GetSolid();
	//G4ThreeVector position1 = G4ThreeVector(0., -4 * mm, 0.);
	//auto solidPENShell = new G4UnionSolid("solidPENShell", solidInnerShell, solidOuterShell, 0, position1);

	auto reflectormesh = CADMesh::TessellatedMesh::FromSTL("../models/OuterShell-II-InnerReflector.stl");
	G4VSolid* solidTemp = reflectormesh->GetSolid();
	auto solidReflector = new G4SubtractionSolid("solidReflector", solidTemp, solidOuterShell, 0, G4ThreeVector());
	auto logicReflector = new G4LogicalVolume(solidReflector, matLN2, "logicReflector");
	return logicReflector;

}

G4LogicalVolume* PENDetectorConstruction::ConstructInnerShell() {
}

G4LogicalVolume* PENDetectorConstruction::ConstructOuterShell() {

}

G4LogicalVolume* PENDetectorConstruction::ConstructSArBrick() {
	auto SArCrystalmesh = CADMesh::TessellatedMesh::FromSTL("../models/SArCrystal-II.stl");
	auto SArBrickmesh = CADMesh::TessellatedMesh::FromSTL("../models/SArContainer-II-solid.stl");
	G4VSolid* solidSArCrystal = SArCrystalmesh->GetSolid();
	G4VSolid* solidSArBrick = SArBrickmesh->GetSolid();

	auto solidSArContainer = new G4SubtractionSolid("solidSArCrystal", solidSArBrick, solidSArCrystal, 0, G4ThreeVector(0, 1 * cm, 0));
	auto logicSArBrick = new G4LogicalVolume(solidSArBrick, matPEN, "logicSArBrick");
	auto logicSArCrystal = new G4LogicalVolume(solidSArCrystal, matLAr, "logicSArCrystal");
	//auto logicSArContainer = new G4LogicalVolume(solidSArCrystal, matPEN, "logicSArContainer");

	auto physSArCrystal = new G4PVPlacement(0, G4ThreeVector(0, 1 * cm, 0), logicSArCrystal, "SArCrystal", logicSArBrick, false, 0, CheckOverlaps);
	//auto physSArContainer = new G4PVPlacement(0, G4ThreeVector(), logicSArContainer, "SArContainer", logicSArBrick, false, 0, CheckOverlaps);
	return logicSArBrick;
}

G4LogicalVolume* PENDetectorConstruction::ConstructSiPMArray() {
	G4double SiPMThickness = 0.5 * mm;
	G4Box* solidSiPMArray = new G4Box("solidSiPMArray", 3 * cm, 3 * cm, SiPMThickness);
	G4LogicalVolume* logicSiPMArray = new G4LogicalVolume(solidSiPMArray, matLN2, "logicSiPMArray");

	G4double SiPMWidth = 0.8 * cm;
	G4double SiPMLength = 1 * cm;
	G4Box* solidSiPM = new G4Box("solidSiPMArray", SiPMWidth, SiPMLength, SiPMThickness);
	G4LogicalVolume* logicSiPM = new G4LogicalVolume(solidSiPM, matSi, "logicSiPM");
	G4double offset1 = SiPMWidth + 0.5 * mm;
	G4double offset2 = SiPMLength + 0.5 * mm;

	physSiPM0 = new G4PVPlacement(0, G4ThreeVector(offset1, offset2, 0), logicSiPM, "physSiPM0", logicSiPMArray, false, 0, CheckOverlaps);
	physSiPM1 = new G4PVPlacement(0, G4ThreeVector(-offset1, offset2, 0), logicSiPM, "physSiPM1", logicSiPMArray, false, 1, CheckOverlaps);
	physSiPM2 = new G4PVPlacement(0, G4ThreeVector(offset1, -offset2, 0), logicSiPM, "physSiPM2", logicSiPMArray, false, 2, CheckOverlaps);
	physSiPM3 = new G4PVPlacement(0, G4ThreeVector(-offset1, -offset2, 0), logicSiPM, "physSiPM3", logicSiPMArray, false, 3, CheckOverlaps);
	return logicSiPMArray;
}

G4LogicalVolume* PENDetectorConstruction::ConstructSArSiPM() {

}

G4VPhysicalVolume* PENDetectorConstruction::ConstructUnit()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  //Vaccum for world
  //G4Material* vacuum=new G4Material("Galactic",z=1.,a=1.01*g/mole,density=universe_mean_density,kStateGas,2.73*kelvin,3.e-18*pascal);

  //------------------------------------------------------ volumes
  //
  G4Material* world_mat = fVacuum;
  G4Material* env_mat = matLN2;
  G4Material* det_mat = matEnGe;

  //     
  // World&Envelope
  //
  G4double world_size = 200 * cm;
  G4double env_size = 180 * cm;
  G4Box* solidWorld = new G4Box("solidWorld", 0.5 * world_size, 0.5 * world_size, 0.5 * world_size);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "logicWorld");
  G4VPhysicalVolume* physWorld =
	  new G4PVPlacement(0,                     //no rotation
		  G4ThreeVector(),       //at (0,0,0)
		  logicWorld,            //its logical volume
		  "World",               //its name
		  0,                     //its mother  volume
		  false,                 //no boolean operation
		  0,                     //copy number
		  checkOverlaps);        //overlaps checking

  G4Box* solidEnv = new G4Box("solidEnvelope",  0.5 * env_size, 0.5 * env_size, 0.5 * env_size);
  G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv, env_mat, "logicEnvelope");
  physEnv = new G4PVPlacement(0, G4ThreeVector(), logicEnv, "Envelope", logicWorld, false, 0, checkOverlaps);

  G4LogicalVolume* logicTotalCrystal = ConstructBEGe();
  auto physDet = new G4PVPlacement(0, G4ThreeVector(), logicTotalCrystal, "PhysDet", logicEnv, false, 0, checkOverlaps);
  
  //========================PEN shell and wire paremeters========================//

  G4String WireType = fWireType;
  G4double wireradius = 0.7 * mm;
  G4double LN2Gap = 0.5 * mm;
  G4double ShellThickness = 1 * cm;
  G4double WireLength = 71 * mm;
  G4double outerGeRadius = 40 * mm;

  //=============================================================================//

  fWireLength = WireLength;
  G4double WireCentDist;
  G4ThreeVector WirePlacement;

  if (WireType == "A1") {
	  G4LogicalVolume* logicWire = ConstructA1(WireLength);
	  fWireCentDist = outerGeRadius + LN2Gap + ShellThickness + fWireRadius;
	  WirePlacement = G4ThreeVector(fWireCentDist, 0, 0);
	  fWirePos = WirePlacement;
	  //physWire = new G4PVPlacement(0, WirePlacement, logicWire, "Wire", logicEnv, false, 0, checkOverlaps);
	  }

  else if (WireType == "A2") {
	  G4LogicalVolume* logicWire = ConstructA2(WireLength);
	  fWireCentDist = outerGeRadius + LN2Gap + ShellThickness + fWireRadius;
	  WirePlacement = G4ThreeVector(fWireCentDist, 0, 0);
	  fWirePos = WirePlacement;
	  //physWire = new G4PVPlacement(0, WirePlacement, logicWire, "Wire", logicEnv, false, 0, checkOverlaps);
	  }

  else {
	  G4cout << "Type does not exist!" << G4endl;
  }

  G4OpticalSurface* Wire_LN2 = new G4OpticalSurface("Wire_LN2");
  G4LogicalBorderSurface* Wire_LN2_LBS = new G4LogicalBorderSurface("Wire_LN2_LBS", physEnv, physWire, Wire_LN2);
  Wire_LN2 = dynamic_cast <G4OpticalSurface*>(Wire_LN2_LBS->GetSurface(physEnv, physWire)->GetSurfaceProperty());
  Wire_LN2->SetType(dielectric_dielectric);
  Wire_LN2->SetModel(DAVIS);
  Wire_LN2->SetFinish(Polished_LUT);
  
  //PEN shell scintilator
  auto rotPENShell = new G4RotationMatrix();
  rotPENShell->rotateX(90 * degree);
  G4LogicalVolume* logicPENShell = ConstructPENShell();
  physPENShell = new G4PVPlacement(rotPENShell, G4ThreeVector(0, 0, 35.1 * mm), logicPENShell, "PENShell", logicEnv, false, 0, checkOverlaps);

  //Reflector
  if (ifOuterReflector == true) {
	  G4LogicalVolume* logicReflector = ConstructOuterReflector();
	  physOuterReflector = new G4PVPlacement(rotPENShell, G4ThreeVector(0, 0, 35.1 * mm), logicReflector, "Reflector", logicEnv, false, 0, checkOverlaps);

	  //Reflector
	  G4OpticalSurface* PEN_OuterReflector = new G4OpticalSurface("PEN_Reflector");
	  G4LogicalBorderSurface* PEN_OuterReflector_LBS = new G4LogicalBorderSurface("PEN_LN2_LBS", physOuterReflector, physPENShell, PEN_OuterReflector);
	  PEN_OuterReflector = dynamic_cast <G4OpticalSurface*>(PEN_OuterReflector_LBS->GetSurface(physOuterReflector, physPENShell)->GetSurfaceProperty());
	  PEN_OuterReflector->SetType(dielectric_LUTDAVIS);
	  PEN_OuterReflector->SetModel(DAVIS);
	  if (fReflectorType == "ESR") {
		  PEN_OuterReflector->SetFinish(PolishedESR_LUT);
	  }
	  else if (fReflectorType == "ESRGrease") {
		  PEN_OuterReflector->SetFinish(PolishedESRGrease_LUT);
	  }
	  else if (fReflectorType == "Teflon") {
		  PEN_OuterReflector->SetFinish(PolishedTeflon_LUT);
	  }
	  else {
		  G4cout << "Reflector type not found! Using default ESR LUT." << G4endl;
		  PEN_OuterReflector->SetFinish(PolishedESR_LUT);
	  }

  }

  if (ifInnerReflector == true) {
	  G4LogicalVolume* logicReflector = ConstructInnerReflector();
	  physInnerReflector = new G4PVPlacement(rotPENShell, G4ThreeVector(0, 0, 35.1 * mm), logicReflector, "Reflector", logicEnv, false, 0, checkOverlaps);

	  //Reflector
	  G4OpticalSurface* PEN_InnerReflector = new G4OpticalSurface("PEN_Reflector");
	  G4LogicalBorderSurface* PEN_InnerReflector_LBS = new G4LogicalBorderSurface("PEN_LN2_LBS", physInnerReflector, physPENShell, PEN_InnerReflector);
	  PEN_InnerReflector = dynamic_cast <G4OpticalSurface*>(PEN_InnerReflector_LBS->GetSurface(physInnerReflector, physPENShell)->GetSurfaceProperty());
	  PEN_InnerReflector->SetType(dielectric_LUTDAVIS);
	  PEN_InnerReflector->SetModel(DAVIS);
	  PEN_InnerReflector->SetFinish(PolishedESR_LUT);
	  if (fReflectorType == "ESR") {
		  PEN_InnerReflector->SetFinish(PolishedESR_LUT);
	  }
	  else if (fReflectorType == "ESRGrease") {
		  PEN_InnerReflector->SetFinish(PolishedESRGrease_LUT);
	  }
	  else if (fReflectorType == "Teflon") {
		  PEN_InnerReflector->SetFinish(PolishedTeflon_LUT);
	  }
	  else {
		  G4cout << "Reflector type not found! Using default ESR LUT." << G4endl;
		  PEN_InnerReflector->SetFinish(PolishedESR_LUT);
	  }
	  //
  }



  //=============================================================//
  //                                                             //
  //                            SiPMs                            //
  //                                                             //
  //=============================================================//

  G4LogicalVolume* logicSiPMArray = ConstructSiPMArray();

  auto rotSiPM0 = new G4RotationMatrix();
  rotSiPM0->rotateY(90 * degree);

  physSiPMArray0 = new G4PVPlacement(rotSiPM0, G4ThreeVector(-167 * mm, 0, 0 * mm), logicSiPMArray, "SiPMArray0", logicEnv, false, 0, checkOverlaps);
  physSiPMArray1 = new G4PVPlacement(rotSiPM0, G4ThreeVector(167 * mm, 0, 0 * mm), logicSiPMArray, "SiPMArray1", logicEnv, false, 1, checkOverlaps);



  G4OpticalSurface* PEN_LN2 = new G4OpticalSurface("PEN_LN2");
  G4OpticalSurface* Ge_LN2 = new G4OpticalSurface("Ge_LN2");
  G4OpticalSurface* SiPM_LN2_0 = new G4OpticalSurface("SiPM_LN2_0");
  G4OpticalSurface* SiPM_LN2_1 = new G4OpticalSurface("SiPM_LN2_1");
  G4OpticalSurface* SiPM_LN2_2 = new G4OpticalSurface("SiPM_LN2_2");
  G4OpticalSurface* SiPM_LN2_3 = new G4OpticalSurface("SiPM_LN2_3");
  G4OpticalSurface* SiPM_LN2_4 = new G4OpticalSurface("SiPM_LN2_4");
  G4OpticalSurface* SiPM_LN2_5 = new G4OpticalSurface("SiPM_LN2_5");

  G4LogicalBorderSurface* PEN_LN2_LBS = new G4LogicalBorderSurface("PEN_LN2_LBS", physEnv, physPENShell, PEN_LN2);
  G4LogicalBorderSurface* Ge_LN2_LBS = new G4LogicalBorderSurface("PEN_LN2_LBS", physEnv, physDet, PEN_LN2);
  G4LogicalBorderSurface* SiPM_LN2_LBS_0 = new G4LogicalBorderSurface("SiPM_LN2_LBS_0", physEnv, physSiPM0, SiPM_LN2_0);
  G4LogicalBorderSurface* SiPM_LN2_LBS_1 = new G4LogicalBorderSurface("SiPM_LN2_LBS_1", physEnv, physSiPM1, SiPM_LN2_1);
  G4LogicalBorderSurface* SiPM_LN2_LBS_2 = new G4LogicalBorderSurface("SiPM_LN2_LBS_2", physEnv, physSiPM2, SiPM_LN2_2);
  G4LogicalBorderSurface* SiPM_LN2_LBS_3 = new G4LogicalBorderSurface("SiPM_LN2_LBS_3", physEnv, physSiPM3, SiPM_LN2_3);
  G4LogicalBorderSurface* SiPM_LN2_LBS_4 = new G4LogicalBorderSurface("SiPM_LN2_LBS_4", physEnv, physSiPM4, SiPM_LN2_4);
  G4LogicalBorderSurface* SiPM_LN2_LBS_5 = new G4LogicalBorderSurface("SiPM_LN2_LBS_5", physEnv, physSiPM5, SiPM_LN2_5);
 // G4LogicalBorderSurface* SiPM_LN2_LBS_P = new G4LogicalBorderSurface("SiPM_LN2_LBS_P", physEnv, physSiPMP, SiPM_LN2_P);

  PEN_LN2 = dynamic_cast <G4OpticalSurface*>(PEN_LN2_LBS->GetSurface(physEnv, physPENShell)->GetSurfaceProperty());
  Ge_LN2 = dynamic_cast <G4OpticalSurface*>(PEN_LN2_LBS->GetSurface(physEnv, physDet)->GetSurfaceProperty());
  SiPM_LN2_0 = dynamic_cast <G4OpticalSurface*>(SiPM_LN2_LBS_0->GetSurface(physEnv, physSiPM0)->GetSurfaceProperty());
  SiPM_LN2_1 = dynamic_cast <G4OpticalSurface*>(SiPM_LN2_LBS_1->GetSurface(physEnv, physSiPM1)->GetSurfaceProperty());
  SiPM_LN2_2 = dynamic_cast <G4OpticalSurface*>(SiPM_LN2_LBS_2->GetSurface(physEnv, physSiPM2)->GetSurfaceProperty());
  SiPM_LN2_3 = dynamic_cast <G4OpticalSurface*>(SiPM_LN2_LBS_3->GetSurface(physEnv, physSiPM3)->GetSurfaceProperty());
  SiPM_LN2_4 = dynamic_cast <G4OpticalSurface*>(SiPM_LN2_LBS_4->GetSurface(physEnv, physSiPM4)->GetSurfaceProperty());
  SiPM_LN2_5 = dynamic_cast <G4OpticalSurface*>(SiPM_LN2_LBS_5->GetSurface(physEnv, physSiPM5)->GetSurfaceProperty());
 // SiPM_LN2_P = dynamic_cast <G4OpticalSurface*>(SiPM_LN2_LBS_P->GetSurface(physEnv, physSiPMP)->GetSurfaceProperty());

  PEN_LN2->SetType(dielectric_LUTDAVIS);
  //PEN_LN2->SetType(dielectric_dielectric);
  Ge_LN2->SetType(dielectric_metal);
  SiPM_LN2_0->SetType(dielectric_LUTDAVIS);
  SiPM_LN2_1->SetType(dielectric_LUTDAVIS);
  SiPM_LN2_2->SetType(dielectric_LUTDAVIS);
  SiPM_LN2_3->SetType(dielectric_LUTDAVIS);
  SiPM_LN2_4->SetType(dielectric_LUTDAVIS);
  SiPM_LN2_5->SetType(dielectric_LUTDAVIS);
 // SiPM_LN2_P->SetType(dielectric_dielectric);

  PEN_LN2->SetModel(DAVIS);
  Ge_LN2->SetModel(glisur);
  SiPM_LN2_0->SetModel(DAVIS);
  SiPM_LN2_1->SetModel(DAVIS);
  SiPM_LN2_2->SetModel(DAVIS);
  SiPM_LN2_3->SetModel(DAVIS);
  SiPM_LN2_4->SetModel(DAVIS);
  SiPM_LN2_5->SetModel(DAVIS);
  //SiPM_LN2_P->SetModel(DAVIS);

  PEN_LN2->SetFinish(Polished_LUT);
  Ge_LN2->SetFinish(polished);
  SiPM_LN2_0->SetFinish(Detector_LUT);
  SiPM_LN2_1->SetFinish(Detector_LUT);
  SiPM_LN2_2->SetFinish(Detector_LUT);
  SiPM_LN2_3->SetFinish(Detector_LUT);
  SiPM_LN2_4->SetFinish(Detector_LUT);
  SiPM_LN2_5->SetFinish(Detector_LUT);
  //SiPM_LN2_P->SetFinish(Detector_LUT);

  const G4int NUMENTRIES_CHIP = 11;
  const double hc = 6.62606876 * 2.99792458 * 100. / 1.602176462;
  G4double sipm_pp[NUMENTRIES_CHIP] = { hc / 600. * eV, hc / 590. * eV, hc / 580. * eV, hc / 570. * eV, hc / 560. * eV, hc / 550. * eV, hc / 540. * eV, hc / 530. * eV, hc / 520. * eV,hc / 510. * eV,hc / 500. * eV };
  G4double sipm_sl[NUMENTRIES_CHIP] = { 0,0,0,0,0,0,0,0,0,0,0 };
  G4double sipm_ss[NUMENTRIES_CHIP] = { 0,0,0,0,0,0,0,0,0,0,0 };
  G4double sipm_bs[NUMENTRIES_CHIP] = { 0,0,0,0,0,0,0,0,0,0,0 };
  G4double sipm_rindex[NUMENTRIES_CHIP] = { 1.406,1.406,1.406,1.406,1.406,1.406,1.406,1.406,1.406,1.406,1.406 };
  G4double sipm_reflectivity[NUMENTRIES_CHIP] = { 0,0,0,0,0,0,0,0,0,0,0 };
  // G4double sipm_efficiency[NUMENTRIES_CHIP] = {0.20,0.21,0.23,0.25,0.26,0.28,0.30,0.32,0.34,0.36,0.38};
  G4double sipm_efficiency[NUMENTRIES_CHIP] = { 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0 };

  G4MaterialPropertiesTable* SIPM_MPT_Surf = new G4MaterialPropertiesTable();
  // SIPM_MPT_Surf->AddProperty("SPECULARLOBECONSTANT",sipm_pp,sipm_sl,NUMENTRIES_CHIP);
  // SIPM_MPT_Surf->AddProperty("SPECULARSPIKECONSTANT",sipm_pp,sipm_ss,NUMENTRIES_CHIP);
  // SIPM_MPT_Surf->AddProperty("BACKSCATTERCONSTANT",sipm_pp,sipm_bs,NUMENTRIES_CHIP);
  SIPM_MPT_Surf->AddProperty("REFLECTIVITY", sipm_pp, sipm_reflectivity, NUMENTRIES_CHIP);
  SIPM_MPT_Surf->AddProperty("EFFICIENCY", sipm_pp, sipm_efficiency, NUMENTRIES_CHIP);
  // SIPM_MPT_Surf->AddProperty("RINDEX",sipm_pp,sipm_rindex,NUMENTRIES_CHIP);

  SiPM_LN2_0->SetMaterialPropertiesTable(SIPM_MPT_Surf);
  SiPM_LN2_1->SetMaterialPropertiesTable(SIPM_MPT_Surf);
  SiPM_LN2_2->SetMaterialPropertiesTable(SIPM_MPT_Surf);
  SiPM_LN2_3->SetMaterialPropertiesTable(SIPM_MPT_Surf);
  SiPM_LN2_4->SetMaterialPropertiesTable(SIPM_MPT_Surf);
  SiPM_LN2_5->SetMaterialPropertiesTable(SIPM_MPT_Surf);

  return physWorld;
}

G4VPhysicalVolume* PENDetectorConstruction::ConstructSArUnit() {

	G4NistManager* nist = G4NistManager::Instance();
	G4bool checkOverlaps = true;

	G4Material* world_mat = fVacuum;
	G4Material* env_mat = matLN2;
	G4Material* det_mat = matEnGe;

	// World&Envelope
	G4double world_size = 450 * cm;
	G4double env_size = 400 * cm;
	G4Box* solidWorld = new G4Box("solidWorld", 0.5 * world_size, 0.5 * world_size, 0.5 * world_size);
	G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "logicWorld");
	G4VPhysicalVolume* physWorld =
		new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(),       //at (0,0,0)
			logicWorld,            //its logical volume
			"World",               //its name
			0,                     //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking

	G4Box* solidEnv = new G4Box("solidEnvelope", 0.5 * env_size, 0.5 * env_size, 0.5 * env_size);
	G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv, env_mat, "logicEnvelope");
	physEnv = new G4PVPlacement(0, G4ThreeVector(), logicEnv, "Envelope", logicWorld, false, 0, checkOverlaps);

	G4LogicalVolume* logicSArBrick = ConstructSArBrick();
	auto rotSArBrick= new G4RotationMatrix();
	rotSArBrick->rotateX(-90 * degree);
	physSArBrick = new G4PVPlacement(rotSArBrick, G4ThreeVector(), logicSArBrick, "SArBrick", logicEnv, false, 0, checkOverlaps);
	
	return physWorld;
}

G4VPhysicalVolume* PENDetectorConstruction::ConstructArray_1() {

}

G4VPhysicalVolume* PENDetectorConstruction::Construct()
{
	if (fPENPropertiesID == 0) {
		fLY = 6000. / MeV;
		absFactor = 1.5;
	}
	else if (fPENPropertiesID == 1) {
		fLY = 3500. / MeV;
		absFactor = 2.58;
	}
	else if (fPENPropertiesID == 2) {
		fLY = 6000. / MeV;
		absFactor = 6.44;
	}
	SetABS(absFactor);
	SetLY(fLY);

	if (fMode == "Unit") {
		return ConstructUnit();
	}
	if (fMode == "SArUnit") {
		return ConstructSArUnit();
	}
	if (fMode == "Array_1") {
		return ConstructArray_1();
	}
	else {
		G4cout << "Error: Mode not fount!" << G4endl;
	}

}
