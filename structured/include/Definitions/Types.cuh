namespace uammd{
namespace structured{

    struct ForceTorque{
        real4 force;
        real4 torque;
    };

    struct EnergyForceTorque{
        real  energy;
        real4 force;
        real4 torque;
    };

    struct ForceTorqueMagneticField{
        real4 force;
        real4 torque;
        real4 magneticField;
    };

    struct StateTransitionProbability{
        int4 tentativeState;
        real transitionProbability;
    };

}}
