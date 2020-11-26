#ifndef PENDetectorConstruction_h
#define PENDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4String.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "PENMaterials.hh"
#include "PENDetectorMessenger.hh"
#include "PENDetectorConstruction.hh"
//#include "TMath.h"
//#include "G4GDMLParser.hh"

class PENMaterials;
class G4VPhysicalVolume;
class G4LogicalVolume;
class PENDeterctorMessenger;

class PENDetectorConstruction : public G4VUserDetectorConstruction
{
    public:
        PENDetectorConstruction();
        ~PENDetectorConstruction();
        const G4VPhysicalVolume* GetPENShell() const;
		const G4VPhysicalVolume* GetBulk() const;
        const G4VPhysicalVolume* GetSiPM(G4int i) const;
        const G4VPhysicalVolume* GetEnv() const;

        G4VPhysicalVolume* Construct();
        G4VPhysicalVolume* ConstructUnit();
        G4VPhysicalVolume* ConstructSArUnit();
        G4VPhysicalVolume* ConstructArray_1();

        G4LogicalVolume* ConstructOuterShell();
        G4LogicalVolume* ConstructInnerShell();
        G4LogicalVolume* ConstructPENShell();
        G4LogicalVolume* ConstructBEGe();
        G4LogicalVolume* ConstructA1(G4double WireLength);
        G4LogicalVolume* ConstructA2(G4double WireLength);
        G4LogicalVolume* ConstructSArBrick();
        G4LogicalVolume* ConstructOuterReflector();
        G4LogicalVolume* ConstructInnerReflector();
        G4LogicalVolume* ConstructSiPMArray();
        G4LogicalVolume* ConstructSArSiPM();

        void DefineMat();

        void SetABS(G4double);
        void SetLY(G4double);
        void SetWireType(G4String);
        void SetReflectorType(G4String);
        void SetConfine(G4String);
        void SetRunInfo(G4String);
        void SetMode(G4String);
        void SetPENPropertiesID(G4int);
        void SetOuterReflector(G4bool);
        void SetInnerReflector(G4bool);

        //void SetLayerNbS(G4String);

        G4ThreeVector GetWirePos() {
            return fWirePos;
        }

        G4double GetWireRadius() {
            return fWireRadius;
        }

        G4double GetWireLength() {
            return fWireLength;
        }

        G4String GetWireType() {
            return fWireType;
        }

        G4String GetConfine() {
            return fConfine;
        }

        G4String GetRunInfo() {
            return fRunInfo;
        }

        G4String GetMode() {
            return fMode;
        }

        G4int GetPENPropertiesID() {
            return fPENPropertiesID;
        }

        G4double GetPENShellLength() {
            return fPENShellLength;
        }

        G4double GetPENShellRadius() {
            return fPENShellRadius;
        }

        G4double GetSArBrickRadius() {
            return fSArBrickRadius;
        }

        G4double GetSArBrickHeight() {
            return fSArBrickHeight;
        }


    private:
        G4VPhysicalVolume* physSArBrick;
		G4VPhysicalVolume* physBulk;
        G4VPhysicalVolume* physEnv;
        G4VPhysicalVolume* physSiPM0;
        G4VPhysicalVolume* physSiPM1;
        G4VPhysicalVolume* physSiPM2;
        G4VPhysicalVolume* physSiPM3;
        G4VPhysicalVolume* physSiPM4;
        G4VPhysicalVolume* physSiPM5;
        G4VPhysicalVolume* physSiPM6;
        G4VPhysicalVolume* physSiPM7;
        G4VPhysicalVolume* physSiPM8;
        G4VPhysicalVolume* physSiPM9;

        G4VPhysicalVolume* physSiPMArray0;
        G4VPhysicalVolume* physSiPMArray1;
        G4VPhysicalVolume* physSiPMArray2;
        G4VPhysicalVolume* physSiPMArray3;

        G4VPhysicalVolume* physWire;
        G4VPhysicalVolume* physPENShell;
        G4VPhysicalVolume* physOuterReflector;
        G4VPhysicalVolume* physInnerReflector;
        //G4LogicalVolume* logicPENShell;
        G4Tubs* solidSideSiPM;

        PENMaterials* matconstructor;

        G4double fLY;
        G4double fRES;
        G4double AbsorptionLength;
        G4double pmtReflectivity;
        G4double fRI;

        G4Material* fWorldMaterial;
        G4Material* fTargetMaterial;
        G4Material* fGlassMaterialPMT;
        G4Material* fPhotoCathodeMaterial;

        G4OpticalSurface* AirTarget;
        G4OpticalSurface* surfaceCathodeSupport;
        G4OpticalSurface* surfaceCathodeSupportBottom;

        G4MaterialPropertiesTable* MPT_PEN;
        G4MaterialPropertiesTable* MPT_GlassPMT;
        G4MaterialPropertiesTable* MPT_Target;
        G4MaterialPropertiesTable* SMPT_AirTarget;
        G4MaterialPropertiesTable* MPT_World;

        G4Material* matPEN;
        G4Material* matBialkali;
        G4Material* matSi;
        G4Material* fGe;
        G4Material* matAir;
        G4Material* fVacuum;
        G4Material* matTriggerFoilEJ212;
        G4Material* Pstyrene;
        G4Material* fGlass;
        G4Material* fPOM;
        G4Material* fABS;
        G4Material* matPMMA;
        G4Material* matEnGe;
        G4Material* matNaGe;
        G4Material* matGreaseEJ550;
        G4Material* matTeflon;
        G4Material* matVikuiti;
        G4Material* matPolyethylene;
        G4Material* matTitanium;
        G4Material* matLN2;
        G4Material* matLAr;
        G4Material* matGAGG;
		G4Material* matPTFE;
		G4Material* matCu;
        G4Material* matNylon;
        G4Material* matTPB;

        G4Material* fDetMat;

        G4String fABSFile;
        G4String fWireType;
        G4String fReflectorType;
        G4String fConfine;
        G4String fRunInfo;
        G4String fMode;
        G4int fPENPropertiesID;
        G4ThreeVector fWirePos;
        G4double fWireCentDist;
        G4double fWireRadius;
        G4double fWireLength;
        G4double fPENShellRadius;
        G4double fPENShellLength;
        G4double absFactor;
        G4double fBEGeRadius;
        G4double fBEGeHeight;
        G4double fSArBrickHeight;
        G4double fSArBrickRadius;
        PENDetectorMessenger* fDetectorMessenger;
        G4bool CheckOverlaps;
        G4bool ifOuterReflector;
        G4bool ifInnerReflector;
        
};

inline const G4VPhysicalVolume* PENDetectorConstruction::GetPENShell() const
{
    return physPENShell;
}

inline const G4VPhysicalVolume* PENDetectorConstruction::GetBulk() const
{
	return physBulk;
}

inline const G4VPhysicalVolume* PENDetectorConstruction::GetEnv() const
{
    return physEnv;
}

inline const G4VPhysicalVolume* PENDetectorConstruction::GetSiPM(G4int i) const
{
    switch  (i){
        case 0:
        return physSiPM0;
        break;
        case 1:
        return physSiPM1;
        break;
        case 2:
        return physSiPM2;
        break;
        case 3:
        return physSiPM3;
        break;
        case 4:
        return physSiPM4;
        break;
        case 5:
        return physSiPM5;
        break;
        case 6:
        return physSiPM6;
        break;
        case 7:
        return physSiPM7;
        break;
        case 8:
        return physSiPM8;
        break;
        case 9:
        return physSiPM9;
        break;
        default:
        break;
    }
}

#endif