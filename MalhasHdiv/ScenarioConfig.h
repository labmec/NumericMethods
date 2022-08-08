#ifndef SCENARIO_CONFIG
#define SCENARIO_CONFIG

enum EScenario {Scenario2, Scenario3, Scenario4, Scenario5, Scenario6, Scenario7, Scenario8, Scenario9};

struct ScenarioConfig
{
    int CoarsePressurePOrder = 1;
    int FinePressurePOrder = 1;
    int CoarseBubbleFluxPOrder = 1;
    int FineBubbleFluxPOrder = 1;
    int CoarseBoundaryFluxPOrder = 1;
    int FineBoundaryFluxPOrder = 1;
    HDivFamily HDivFam;

    //To configurate a scenario, chose its number and set the Base Polynomial Order
    EScenario Scenario;
    int BasePOrder = 1;

    void ConfigurateScenario(){
        switch (Scenario)
        {
        case Scenario2:
            HDivFam = HDivFamily::EHDivStandard;
            CoarsePressurePOrder = 1;
            FinePressurePOrder = BasePOrder;
            CoarseBubbleFluxPOrder = 1;
            FineBubbleFluxPOrder = BasePOrder;
            CoarseBoundaryFluxPOrder = 1;
            FineBoundaryFluxPOrder = BasePOrder;
            break;

        case Scenario3:
            HDivFam = HDivFamily::EHDivStandard;
            CoarsePressurePOrder = 1;
            FinePressurePOrder = 1;
            CoarseBubbleFluxPOrder = 1;
            FineBubbleFluxPOrder = BasePOrder;
            CoarseBoundaryFluxPOrder = 1;
            FineBoundaryFluxPOrder = BasePOrder;
            DebugStop();//This scenario is not consistent
            break;

        case Scenario4:
            HDivFam = HDivFamily::EHDivConstant;
            CoarsePressurePOrder = 0;
            FinePressurePOrder = 0;
            CoarseBubbleFluxPOrder = 1;
            FineBubbleFluxPOrder = BasePOrder;
            CoarseBoundaryFluxPOrder = 1;
            FineBoundaryFluxPOrder = BasePOrder;
            break;

        case Scenario5:
            HDivFam = HDivFamily::EHDivStandard;
            CoarsePressurePOrder = BasePOrder;
            FinePressurePOrder = BasePOrder;
            CoarseBubbleFluxPOrder = BasePOrder;
            FineBubbleFluxPOrder = BasePOrder;
            CoarseBoundaryFluxPOrder = 1;
            FineBoundaryFluxPOrder = BasePOrder;
            break;

        case Scenario6:
            HDivFam = HDivFamily::EHDivStandard;
            CoarsePressurePOrder = 1;
            FinePressurePOrder = 1;
            CoarseBubbleFluxPOrder = BasePOrder;
            FineBubbleFluxPOrder = BasePOrder;
            CoarseBoundaryFluxPOrder = 1;
            FineBoundaryFluxPOrder = BasePOrder;
            DebugStop();//This scenario is not consistent
            break;

        case Scenario7:
            HDivFam = HDivFamily::EHDivConstant;
            CoarsePressurePOrder = 0;
            FinePressurePOrder = 0;
            CoarseBubbleFluxPOrder = BasePOrder;
            FineBubbleFluxPOrder = BasePOrder;
            CoarseBoundaryFluxPOrder = 1;
            FineBoundaryFluxPOrder = BasePOrder;
            break;
        
        case Scenario8:
            HDivFam = HDivFamily::EHDivConstant;
            CoarsePressurePOrder = 0;
            FinePressurePOrder = 0;
            CoarseBubbleFluxPOrder = 1;
            FineBubbleFluxPOrder = BasePOrder;
            CoarseBoundaryFluxPOrder = 0;
            FineBoundaryFluxPOrder = BasePOrder;
            break;

        case Scenario9:
            HDivFam = HDivFamily::EHDivConstant;
            CoarsePressurePOrder = 0;
            FinePressurePOrder = 0;
            CoarseBubbleFluxPOrder = BasePOrder;
            FineBubbleFluxPOrder = BasePOrder;
            CoarseBoundaryFluxPOrder = 0;
            FineBoundaryFluxPOrder = BasePOrder;
            break;

        default:
            DebugStop();
            break;
        }
    }
};



#endif
