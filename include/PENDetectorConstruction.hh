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
#include <map>
//#include "TMath.h"
//#include "G4GDMLParser.hh"

class PENMaterials;
class G4VPhysicalVolume;
class G4LogicalVolume;
class PENDeterctorMessenger;
class G4AssemblyVolume;

class PENDetectorConstruction : public G4VUserDetectorConstruction
{
    public:
        PENDetectorConstruction();
        ~PENDetectorConstruction();
        const G4VPhysicalVolume* GetPENShell() const;
		const G4VPhysicalVolume* GetBulk() const;
        const G4VPhysicalVolume* GetSiPM(G4int i) const;
        const G4VPhysicalVolume* GetContainerSiPM(G4int i) const;
        const G4VPhysicalVolume* GetEnv() const;

        G4VPhysicalVolume* Construct();
        G4VPhysicalVolume* ConstructUnit();
        G4VPhysicalVolume* ConstructSArUnit();
        G4VPhysicalVolume* ConstructArray_1();
        G4VPhysicalVolume* GetPhysicalVolumeByName(const G4String& name);
        void ResetPhysicalVolumeNames();
        void GetPhysicalVolumeProperties();
        std::map<G4VPhysicalVolume*, G4int> VolumeLUT;
        std::map<G4VPhysicalVolume*, G4String> VolumeNameLUT;

        G4LogicalVolume* ConstructOuterShell();
        G4LogicalVolume* ConstructInnerShell();
        G4LogicalVolume* ConstructPENShell();
        G4LogicalVolume* ConstructCSGOuterShell();
        G4LogicalVolume* ConstructCSGInnerShell();
        G4LogicalVolume* ConstructCSGPENShell();
        G4LogicalVolume* ConstructBEGe();
        G4LogicalVolume* ConstructA1(G4double WireLength);
        G4LogicalVolume* ConstructA2(G4double WireLength);
        G4LogicalVolume* ConstructContainerBrick();
        G4LogicalVolume* ConstructStringBoxBrick();
        G4LogicalVolume* ConstructStringBox();
        G4LogicalVolume* ConstructOuterReflector();//Outer Reflector of OuterShell
        G4LogicalVolume* ConstructInnerReflector();//Inner Reflector of OuterShell
        G4LogicalVolume* ConstructReflector();
        G4LogicalVolume* ConstructASICPlate();

        G4LogicalVolume* ConstructSArSiPMArrayLV();
        G4LogicalVolume* ConstructContainerSiPMArrayLV();
        G4LogicalVolume* ConstructSiPMArrayLV();

        G4AssemblyVolume* ConstructSArSiPMArray();
        G4AssemblyVolume* ConstructContainerSiPMArray();
        G4AssemblyVolume* ConstructSiPMArray();


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
        G4String GetMode() {
            return fMode;
        }

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

        G4String GetReflectorType() {
            return fReflectorType;
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
        G4double GetASICThickness() {
            return fASICThickness;
        }


    private:
        G4VPhysicalVolume* physContainerBrick;
        G4VPhysicalVolume* physStringBoxBrick;//String Box Brick
		G4VPhysicalVolume* physBulk;
        G4VPhysicalVolume* physEnv;
        G4VPhysicalVolume* physSiPM0;
        G4VPhysicalVolume* physSiPM1;
        G4VPhysicalVolume* physSiPM2;
        G4VPhysicalVolume* physSiPM3;
        G4VPhysicalVolume* physSiPM4;
        G4VPhysicalVolume* physSiPM5;

        G4VPhysicalVolume* physContainerSiPM0;
        G4VPhysicalVolume* physContainerSiPM1;
        G4VPhysicalVolume* physContainerSiPM2;
        G4VPhysicalVolume* physContainerSiPM3;
        G4VPhysicalVolume* physContainerSiPM4;
        G4VPhysicalVolume* physContainerSiPM5;

        G4VPhysicalVolume* physSiPMArray0;
        G4VPhysicalVolume* physSiPMArray1;
        G4VPhysicalVolume* physSiPMArray2;
        G4VPhysicalVolume* physSiPMArray3;

        G4VPhysicalVolume* physWire;
        G4VPhysicalVolume* physPENShell;
        G4VPhysicalVolume* physInnerShell;
        G4VPhysicalVolume* physOuterShell;
        G4VPhysicalVolume* physOuterReflector;
        G4VPhysicalVolume* physInnerReflector;
        G4VPhysicalVolume* physContainerCrystal;//Crystal in Container
        G4VPhysicalVolume* physStringBoxCrystal;//Crystal in String Box
        G4VPhysicalVolume* physASICPlate;
        G4LogicalVolume* logicStringBoxCrystal;
        G4LogicalVolume* logicContainerCrystal;
        G4LogicalVolume* logicASICPlate;
        G4LogicalVolume* logicWire;

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
        G4double fASICWidth;
        G4double fASICLength;
        G4double fASICThickness;
        G4double fOuterShellLength;
        G4double fOuterShellHeight;
        G4double fOuterShellWidth;

        PENDetectorMessenger* fDetectorMessenger;
        G4bool CheckOverlaps;
        G4bool ifOuterReflector;
        G4bool ifInnerReflector;
        G4bool ifReflector;
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
        default:
        break;
    }
}

inline const G4VPhysicalVolume* PENDetectorConstruction::GetContainerSiPM(G4int i) const
{
    switch (i) {
    case 0:
        return physContainerSiPM0;
        break;
    case 1:
        return physContainerSiPM1;
        break;
    case 2:
        return physContainerSiPM2;
        break;
    case 3:
        return physContainerSiPM3;
        break;
    case 4:
        return physContainerSiPM4;
        break;
    default:
        break;
    }
}

#endif