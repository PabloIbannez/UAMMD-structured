Design
======
Simulation
----------

UAMMD-structured compiles into a single executable. However, it can also
be used as a library by utilizing specific parts of the code or
the entirety of it. To include the entire UAMMD-structured, we can use
``#include "UAMMDstructured.cuh"``. This is precisely what we do when
generating the main binary:

.. code-block:: cpp

   #include "UAMMDstructured.cuh"

    using namespace uammd::structured;

    int main(int argc, char *argv[]) {
        //Check if the input parameters are correct
        if (argc < 2) {
            uammd::System::log<uammd::System::CRITICAL>("No input file provided!");
            return EXIT_FAILURE;
        }
        //Init system. System takes the input file path as argument
        std::string inputFilePath = argv[1];
        std::shared_ptr<ExtendedSystem> sys = std::make_shared<ExtendedSystem>(inputFilePath);
        {
            //Create simulation object
            std::shared_ptr<Simulation> sim = std::make_shared<Simulation>(sys);
            //Run the main simulation loop
            sim->run();
        }
        //End simulation
        sys->finish();
        return EXIT_SUCCESS;
    }

The code above is the way the simulation is executed. First, we
initialize *System*, which is then passed as an argument to the
*Simulation* class, which we will describe next. As observed, *System*
takes the path to the input file as an argument. System will parse the
file, and if correct, the rest of the classes will be able to request it
via ``sys->getInput()`` and access its information.

The initialization process is as follows:

.. code-block:: cpp

    class Simulation{
    private:
       ...
    public:
       ...
       Simulation(std::shared_ptr<ExtendedSystem> sys):sys(sys){

           // Load global data
           gd  = std::make_shared<GlobalData>(sys);
           // Load particle data
           pd  = std::make_shared<ExtendedParticleData>(sys);

           //Load topology
           topology = std::make_shared<Topology>(gd, pd);

           //Load force field
           ff = std::make_shared<ForceField>(topology);

           //Load integrators
           //Note that the integrators are loaded after the topology and the force field
           //this is because topology can set some particle properties that are needed
           //by integrators initialization.
           //For example, the particle mass and the particle radius.
           integrators = std::make_shared<IntegratorManager>(topology);

           //Load simulations steps
           simulationSteps = std::make_shared<SimulationStepManager>(integrators,ff);

           //Handle backup
           ...
       }
       ...
    };

Firstly, we load the
structures that store the data, namely, ``GlobalData`` and
``ParticleData``\  [1]_. GlobalData processes the information in the
*global* section, and the *state* section will be processed by
``ParticleData``.

Subsequently, we use global data and particle data to load the topology.
In the topology, we process the *structure* section, where particle
types and the structures they belong to, like the molecule, residue,
etc., are defined. Also, in topology, we process the interactions
present in the forcefield section. When we process forcefield, *state*,
*global*, and *structure* have already been processed, so the
interactions can already assume this. We then extract the forcefield
object from the topology. In forcefield, we process the interactors to
get them all working together. The forcefield is itself an interactor.
In fact, it is an interactor that aggregates a set of interactors. Next,
we use topology to initialize the ``IntegratorManager``. Finally, we use
*integrators* and *topology* to initialize the simulation steps. The
reason that the simulation steps need these two objects is that it
endows them with the capability to calculate energies, forces, etc.,
which is often necessary to calculate physical quantities of interest.
This concludes the initialization process of the simulation [2]_

Once all the necessary elements have been initialized, we can start the
simulation. This is done via the ``run`` method of the ``Simulation``
class:

.. code-block:: cpp

    class Simulation{
    private:
       ...
    public:
       ...
       int run(){
          Timer tim; // Timer to measure the time taken by the simulation
          tim.tic(); // Start the timer
          std::map<std::string,
                   std::shared_ptr<SimulationStep::SimulationStepBase>> 
                   simSteps = simulationSteps->getSimulationSteps();
          // Get the simulation steps are requested from the SimulationStepManager 
          // and stored in a map
          System::log<System::MESSAGE>("[Simulation] Running simulation...");
          for(auto& integratorInfo: integrators->getSortedIntegratorSteps()){
              // Iterate through the integrators
              // Get the integrator name and the number of steps for the integrator
              std::string name  = integratorInfo.name;
              ullint      steps = integratorInfo.steps;
              System::log<System::MESSAGE>("[Simulation] Running integrator (%s)"
                                           " for %llu steps...", name.c_str(), 
                                           steps);
              // Get the integrator from the IntegratorManager
              std::shared_ptr<Integrator> currentIntegrator = 
              integrators->getIntegrator(name);
              // Load the force field into the integrators
              currentIntegrator->addInteractor(ff);
              // Initialize the simulation steps
              System::log<System::DEBUG>("[Simulation] Initializing simulation "
                                         "steps...");
              for(auto sStep : simSteps){sStep.second->tryInit();}
              // Iterate through the steps
              for(ullint i = 0; i < steps; i++){
                  for(auto sStep : simSteps){
                      // At each step, apply the simulation steps
                      sStep.second->tryApplyStep();
                  }
                  // Move the integrator forward in time e.i. integrate
                  currentIntegrator->forwardTime();
              }
          }
          // Stop the timer and get the total time taken by the simulation
          auto totalTime = tim.toc();
          // Print the mean time per step
          real fps = real(gd->getFundamental()->getCurrentStep())/totalTime);
          System::log<System::MESSAGE>("[Simulation] Mean FPS: %f",fps);
          return 0; // Successful completion
     }
    }

The run begins by initializing a timer, which will tell us at the end of
the simulation how much time has passed and will be used to calculate
the number of steps (frames) per second (FPS). We then store in a map
all the simulation steps present in the simulation and then begin to
iterate over each of the integrators present. For each integrator, we
will perform a series of steps. Before starting the simulation itself,
we load the forcefield (which, as mentioned, is an interactor) into the
current integrator. This is done via the standard UAMMD interface
``addInteractor(std::shared_ptr<uammd::Interactor> interactor)``. Once
the forcefield has been loaded, we prepare the simulation steps and
start the simulation. In each step, we try to apply the simulation
steps; the simulation step itself accesses the current step via global
data and applies itself according to its internal rules. Once all the
integration steps of all the integrators are completed, we stop the
timer and calculate the FPS, concluding the simulation.




System
------

In UAMMD-structured, the initial object to be initialized, and one that
is unique in its instance, is the *System*. This component acts as a
versatile container, encapsulating functionalities that do not
comfortably fit within other objects. Due to its role as the
first-created object, *System* is readily accessible in many parts of
the code, making it a central element in UAMMD-structured’s
architecture.

*System* stands as a universally available object, as it is a
prerequisite for other objects, which typically take it as an argument.
Furthermore, these objects offer a method, ``getSystem()``, for easy
access. In UAMMD-structured, *System* not only carries out the
responsibilities it has in the standard UAMMD framework but also manages
two critical areas: the handling of the simulation’s input and
overseeing the state and backup processes.

The extended functionality of *System* in UAMMD-structured is
encapsulated in the class ``ExtendedSystem_``, which inherits from
UAMMD’s `System class <https://uammd.readthedocs.io/en/latest/System.html>`_. The class is defined with a template parameter
``InputType_`` and introduces an enumeration ``SIMULATION_STATE`` to
signify the running or stopped state of the simulation.

.. code-block:: cpp

    template<class InputType_>
    class ExtendedSystem_ : public uammd::System {
        public:
            using InputType = InputType_;
            //Enum available states, running and stopped
            enum SIMULATION_STATE {RUNNING, STOPPED};
        private:
            //Simulation state
            SIMULATION_STATE state = RUNNING;
            //Smart pointer to input
            std::shared_ptr<InputType> input;
            //Other attributes (name,backup...)
            ...
            void loadSimulationInformation(std::string entryName){...}
            void loadSimulationBackup(std::string entryName){...}
            void init(){
                ...
                //Iterate over entries in system section
                for(std::string entryName : input->getEntriesList(this->path)){
                    //If entry is Information, loadSimulatonInformation(entryName)
                    //If entry is Backup, loadSimulatonBackup(entryName)
                    //Else, emit an error.
                    ...
                }
                //Check if Information was loaded, else emit an error.
            }
        public:
            ExtendedSystem_(int argc, char *argv[],
                            std::string inputFilePath,
                            std::vector<std::string> path):
                            //Path refer to path in input file (JSON)
                            //By default, path = "system"
                            uammd::System(argc,argv),
                            path(path){
                //Check if the input file exists
                ...
                try{
                    input=std::make_shared<InputType>(inputFilePath);
                }catch(std::exception &e){
                    //Emit error, bad input file
                }
                this->init();
            }
            std::shared_ptr<InputType> getInput(){return input;}
            //Attributes getters and setters (STATE, seed, backup ...)
            ...
    };

    using ExtendedSystem = ExtendedSystem_<InputJSON::InputJSON>;

``ExtendedSystem_`` is customized to accommodate UAMMD-structured’s
specific requirements, including input processing and backup management.
This extended version of *System* ensures that the class is not only
backward compatible with UAMMD but also aligns with the additional
functionalities unique to UAMMD-structured. The use of polymorphism in
C++ allows ``ExtendedSystem_`` to be used interchangeably with UAMMD’s
*System*.

*System* is integral to UAMMD-structured’s operation, responsible for
monitoring the simulation’s state, which oscillates between RUNNING and
STOPPED. A change to the STOPPED state triggers the cessation of the
simulation. This feature is instrumental in terminating the simulation
under certain predefined conditions. Given the system’s ubiquitous
access across the framework, any object can intervene to cease the
simulation by adjusting the *System*\ ’s state.

Apart from state management, *System* also oversees the backup
variables, facilitating a crucial aspect of simulation resilience.

Global
------

‘GlobalData‘ is essentially a class acting as a container for other
class instances managing ‘Fundamental‘, ‘Units‘, ‘Types‘, and
‘Ensemble‘. For its initialization, ‘GlobalData‘ requires only a
‘System‘ instance, which it transmits to the other handler classes,
mainly to access the input and read the necessary parameters.

.. code-block:: cpp

    class GlobalData{
       private:
            std::shared_ptr<ExtendedSystem> sys;
            //Handlers objects
            std::shared_ptr<Units::UnitsHandler>             unitsHandler;
            std::shared_ptr<Fundamental::FundamentalHandler> fundamentalHandler;
            std::shared_ptr<Types::TypesHandler>             typesHandler;
            std::shared_ptr<Ensemble::EnsembleHandler>       ensembleHandler;
            void init(){
                //With InputEntryMaknager we can access all the entries
                //in the "path" section of input. In this case path = "global"
                globalInfo = std::make_shared<InputEntryManager>(sys,path);
                //Try to load fundamental
                if(globalInfo->isEntryPresent("fundamental")){
                    fundamentalHandler =
                    FundamentalLoader::loadFundamental(sys,fundamentalPath);
                }else{
                    //If entry fundamental is not present,
                    //we add a fundamental entry of type "Time"
                    //to input.
                    //Create fundamentalHandler
                }
                //Try to load units
                if(globalInfo->isEntryPresent("units")){
                    unitsHandler =
                    UnitsLoader::loadUnits(sys,unitsPath);
                }else{
                    //If entry units is not present,
                    //we add a units entry of type "None"
                    //to input.
                    //Create unitsHandler
                }
                //Load types and ensemble.
                //If types or ensemble are not present
                //error is emited and simulation ends
            }
        public:
            GlobalData(std::shared_ptr<ExtendedSystem>  sys,
                       //Path refer to path in input file (JSON)
                       //By default, path = "global"
                       std::vector<std::string>         path):path(path){
                this->init();
            }
            std::shared_ptr<ExtendedSystem> getSystem(){...}

            std::shared_ptr<Units::UnitsHandler>             getUnits()      {...}
            std::shared_ptr<Ensemble::EnsembleHandler>       getEnsemble()   {...}
            std::shared_ptr<Types::TypesHandler>             getTypes()      {...}
            std::shared_ptr<Fundamental::FundamentalHandler> getFundamental(){...}
    };

Once the various handlers are initialized, they provide access to the
stored information. Let’s examine each of these handlers separately.
UAMMD-structured incorporates a unit manager. Setting a unit system
essentially amounts to fixing the value of certain constants. This is
implemented by making ‘UnitsHandler‘ a virtual class with public methods
as follows:

.. code-block:: cpp

    class UnitsHandler{
        protected:
            //SubType will store the units name
            std::string subType;
        public:
            UnitsHandler(DataEntry& data){
                subType = data.getSubType();
            }
            ...
            virtual real getBoltzmannConstant(){
                System::log<System::CRITICAL>(
                "[Units] BoltzmannConstant not defined for units \"%s\".",
                subType.c_str());
            }
            virtual real getElectricConversionFactor(){
                System::log<System::CRITICAL>("[Units] ElectricConversionFactor not defined for units \"%s\".",
                subType.c_str());
            }
            ...
    };

.. code-block:: cpp

    class KcalMol_A: public UnitsHandler{
        public:
            KcalMol_A(DataEntry& data):UnitsHandler(data){}

            real getBoltzmannConstant()        override {return 1.987191E-03;}
            real getElectricConversionFactor() override {return 332.0716;}
    };

When setting a unit system, we derive from this class and override the
methods to return the values of these constants for the particular
system. For instance, in the :math:`(\text{Kcal/mol})/\text{Å}` unit
system, the electric conversion factor would be :math:`332.0716`, and
Boltzmann’s constant would be :math:`1.987191E-03`.

Using polymorphism, when we assign a specific units object to
‘unitsHandler‘ (such as "KcalMol_A"), methods like
``getBoltzmannConstant()`` will return values pertinent to the specific
system.

In the ‘GlobalData‘ code, the assignment of ‘unitsHandler‘ is done using
the auxiliary function ``loadUnits``, which takes a ‘System‘ instance
and "unitsPath" as arguments. "unitsPath" is a vector of std::string
containing the path to the units information in the input; by default
"unitsPath" = ["global","units"]. The function is as follows:

.. code-block:: cpp

    std::shared_ptr<typename Units::UnitsHandler>
    loadUnits(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path){

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string unitsType    = data.getType();
        std::string unitsSubType = data.getSubType();

        std::shared_ptr<typename Units::UnitsHandler> units;
        bool found = false;
        if("Units" == unitsType and "None" == unitsSubType){
            System::log<System::MESSAGE>(
            "[UnitsLoader] (%s) Detected None units",
            path.back().c_str());
            units = std::make_shared<Units::None>(data);
            found = true;
        }
        if("Units" == unitsType and "KcalMol_A" == unitsSubType){
            System::log<System::MESSAGE>(
            "[UnitsLoader] (%s) Detected KcalMol_A units",
            path.back().c_str());
            units = std::make_shared<Units::KcalMol_A>(data);
            found = true;
        }
        if(not found){
            System::log<System::CRITICAL>(
            "[UnitsLoader] (%s) Could not find units %s::%s",
            path.back().c_str(),unitsType.c_str(),unitsSubType.c_str());
        }
        return units;
    }

‘Ensemble‘ and ‘Fundamental‘ work similarly, defining certain functions
that are overwritten based on the type of ensemble or fundamental. For
example, in ‘Ensemble‘, we define the ‘getBox‘ function, which returns a
box. This function is overwritten for the NVT ensemble, but would emit
an error in another ensemble where it does not apply.

The handling of ‘Types‘ is slightly different. Initially, we define a
base class ‘TypesHandler‘, which processes the input and stores the
information in a set of dictionaries. Internally, each type is
associated with an integer. ‘TypesHandler‘ has a virtual method that
takes a ‘ParticleData‘ instance and applies the information for each
type.

An auxiliary class ‘Types\_‘ is declared, inheriting from
‘TypesHandler‘. This class requires a template argument for the specific
type, as shown in the example with ‘Basic‘. Each type must have two
static methods: ‘loadType‘ and ‘loadTypesIntoParticleData‘. ‘loadType‘
processes the input and populates the ‘nameToData‘ dictionary, taking
‘typeData‘ as an argument. The function ‘loadTypesIntoParticleData‘
iterates over the particles in ‘ParticleData‘, associating relevant
variables based on their type. For example, ‘Basic‘ associates mass,
radius, and charge to each type. Thus, it iterates over each particle,
setting these variables for each one.

.. code-block:: cpp

    class TypesHandler{
        protected:
            std::map<int,std::string> idToName;
            std::map<std::string,int> nameToId;
            std::map<std::string,std::map<std::string,real>> nameToData;
        public:
            TypesHandler(DataEntry& data){
                auto typesData = data.getDataMap();
                //Load type data, set up idToName,nomeToId adn nameToData
            }
            ...
            virtual void
            loadTypesIntoParticleData(std::shared_ptr<ParticleData> pd) = 0;
    };
    template<class T>
    class Types_: public TypesHandler{
        public:
            Types_(DataEntry& data): TypesHandler(data){
                auto typesData = data.getDataMap();
                for(auto& type: typesData){
                    try{T::loadType(this->nameToData,type);}
                    catch (std::exception& e){/*Emit error*/}
                }
            }
            void loadTypesIntoParticleData(std::shared_ptr<ParticleData> pd)
            override {
                T::loadTypesIntoParticleData(pd,this->idToName,this->nameToData);
            }
    };
    //Example from Basic.cuh:
    struct Basic_{
        template<typename T>
        static void loadType(
                    std::map<std::string,std::map<std::string,real>>& nameToData,
                    std::map<std::string,T>& typeData){

            std::string name = typeData.at("name");
            nameToData[name]["mass"]   = real(typeData.at("mass"));
            ...
        }
        static void loadTypesIntoParticleData(
                    std::shared_ptr<ParticleData> pd,
                    std::map<int,std::string>&    idToName,
                    std::map<std::string,std::map<std::string,real>>& nameToData){
            //Iterate over all particles and load type information (mass,radius and charge for Basic)
        }
    };
    using Basic = Types_<Basic_>;

Note that when ‘Global‘ is declared in the system, particles have not
been added or associated with types yet; this is done subsequently. It
is at this point that ‘Global‘ requires the types, and the function
‘loadTypesIntoParticleData‘ is called.

State
-----

In a manner similar to handling the ‘System‘, we work with an object
derived from UAMMD, specifically an object derived from `ParticleData <https://uammd.readthedocs.io/en/latest/ParticleData.html>`_,
which we call ‘ExtendedParticleData‘. ‘ExtendedParticleData‘ is similar
to ‘ParticleData‘, but it initializes by accessing the input to process
and load the content of the ‘State‘ Data Entry into ‘ParticleData‘:

.. code-block:: cpp

    class ExtendedParticleData: public uammd::ParticleData{
        private:
            std::vector<std::string> path;
        public:

            ExtendedParticleData(std::shared_ptr<ExtendedSystem> sys,
                                 std::vector<std::string>       path):
            uammd::ParticleData(sys->getInput()->getDataEntry(path).getDataSize(),
                                sys),
            path(path){
                auto data = this->getSystem()->getInput()->getDataEntry(path);
                stateLoader(this, data);
            }
            ...
    };

The function ‘StateLoader‘ is responsible for loading the state
information. It does this by checking all available options, and if an
option is present in the labels, it requests this information from the
input and loads it into each particle:

.. code-block:: cpp

    void stateLoader(ParticleData* pd,DataEntry& data){
        std::vector<std::string>   labels = data.getLabels();
        ...
        //Id and position are compulsory
        //Check id label is present
        if(std::find(labels.begin(), labels.end(), "id") == labels.end()){
            System::log<System::CRITICAL>(
            "[StateLoader] Label 'id' not found in the state file."
            );
        } else {
            //Load ids
            ...
        }
        //Check position label is present
        if(std::find(labels.begin(), labels.end(), "position") == labels.end()){
            System::log<System::CRITICAL>(
            "[StateLoader] Label 'position' not found in the state file."
            );
        } else {
            //Load position
            ...
        }
        if(std::find(labels.begin(), labels.end(), "velocity") != labels.end()){
            //Load velocity
        }
        if(std::find(labels.begin(), labels.end(), "direction") != labels.end()){
            //Load direction
        }
        ...
        //Load other avaible state variables
        ...
        //Check all labels have been loaded properly
    }

Two fields are always required: "id" and "position", which respectively
indicate the unique identifier for each particle during the simulation
(a particle’s id remains constant throughout) and its position.

Expanding the capabilities of ‘StateLoader‘ is relatively
straightforward. This involves adding a new label to the list of
available labels and loading its content if it is present in the input.

Integrators
-----------

While `Integrator <https://uammd.readthedocs.io/en/latest/Integrator/index.html>`_ is a standard UAMMD class, UAMMD-structured
necessitates a class to manage various integrators and handle the
schedule. This role is fulfilled by ‘IntegratorManager‘. Upon
initialization, this class invokes the functions ``void loadSchedule()``
and ``void loadIntegrators()``. The former function searches for and
processes the ‘schedule‘ entry, populating the ‘integratorSteps‘
dictionary, an attribute of ‘IntegratorManager‘. Subsequently, the
method ``void loadIntegrators()`` is called. This method iterates over
the ‘integrators‘ section, attempting to load all entries except
‘schedule‘. Loading is facilitated by
``IntegratorLoader::loadIntegrators(...)``, which returns a pointer to the
integrator interface from which all integrators are derived. This
process links the type and subtype pair to an instance of an integrator
(a class derived from ‘Integrator‘). Integrators must be predefined in
the schedule and have a unique name; otherwise, an error is thrown,
terminating the program. These integrators are stored in the
‘integrators‘ dictionary.


.. code-block:: cpp

    class IntegratorManager{
      private:
         ...
         std::map<std::string,std::shared_ptr<Integrator>> integrators;
         struct stepsInfo{
             std::string name;
             uint order;
             ullint steps;
         };
         std::map<std::string,stepsInfo> integratorSteps;
         ...
         void loadSchedule(){ /*Set up integratosSteps from schedule entry*/ }
         ...
         void loadIntegrators(){
           for(auto& entry : integratorsInfo->getEntriesInfo()){
             if(IntegratorLoader::isIntegratorAvailable(sys,entry.second.path)){
                 std::shared_ptr<Integrator> integrator = 
                 IntegratorLoader::loadIntegrator(sys,gd,
                                                  groups,
                                                  entry.second.path);
                 if(integrators.count(entry.second.name) == 0){
                     //Check if integrator is in schedule
                     integrators[entry.second.name] = integrator;
                 }
                 else{
                     //Emit error
                 }
             }
           }
           ...
         }
      public:
        IntegratorManager(std::shared_ptr<ExtendedSystem>       sys,
                          std::shared_ptr<GlobalData>            gd,
                          std::shared_ptr<ExtendedParticleData>  pd,
                          std::vector<std::string> path):
                          sys(sys),gd(gd),pd(pd),path(path){
            ...
            this->loadSchedule();
            this->loadIntegrators();
            ...
        }
        ...
        std::vector<stepsInfo> getSortedIntegratorSteps(){...}
        ...
        std::shared_ptr<Integrator> getIntegrator(std::string name){
            ...
            return integrators[name];
        }
    };

‘IntegratorManager‘ possesses several methods to access stored
information. ``getSortedIntegratorSteps()`` returns a list of integrator
information (name, order, and number of steps) sorted by order.
``std::shared_ptr<Integrator> getIntegrator(std::string name)``
retrieves an integrator instance by its name. The integration process
then follows: the list of integrators is retrieved, iterated upon, the
relevant integrator is requested from ‘IntegratorHandler‘, and then
integration occurs for a specified number of steps.

Topology
--------

.. code-block:: cpp

    class Topology{
      private:
        ...
        std::map<std::string,
                 std::shared_ptr<VerletConditionalListSetBase>> VConListSet;
        std::map<std::string,
                 std::shared_ptr<typename uammd::Interactor>> interactors;
        void loadStructure(){...}
        // Methods for force field section processing
        ...
        void loadNeighbourLists(){...}
        void loadInteractors(){
          for(auto& entry : forceFieldInfo->getEntriesInfo()){
           if(Potentials::GenericLoader::isInteractorAvailable(sys,
                                                               entry.second.path)){
              std::shared_ptr<typename uammd::Interactor> inter = 
              Potentials::GenericLoader::loadGeneric(sys,
                                                     gd,groups,
                                                     VConListSet,
                                                     entry.second.path);
              if(interactors.count(entry.second.name) == 0){
                  interactors[entry.second.name] = inter;
                  ...
              } else {/*Emit error*/}
           }
          }
          ...
        }
      public:
        Topology(std::shared_ptr<ExtendedSystem>       sys,
                 std::shared_ptr<GlobalData>            gd,
                 std::shared_ptr<ExtendedParticleData>  pd,
                 std::vector<std::string> path):sys(sys),gd(gd),pd(pd),path(path){
            this->loadStructure();
            //Load components
            ...
            this->loadNeighbourLists();
            this->loadInteractors();
        }
        //Add a new interactor to the system
        void addInteractor(std::shared_ptr<typename uammd::Interactor> interactor, 
                           std::string name){
            ...
        }
        //Getters
        //Get neighbout list std::shared_ptr<VerletConditionalListSetBase>
        //Get interactor std::shared_ptr<typename uammd::Interactor>, by name, type ...
        ...
    };

Structure
^^^^^^^^^

The structure’s loading process is quite similar to that of the ‘State‘.
It involves checking the available labels to see which ones are active
and then loading them accordingly. The process for types, however, is
slightly different.

.. code-block:: cpp

    class Topology{
      private:
        ...
        void loadStructure(){
          ...
          auto structureData = sys->
                               getInput()->
                               getDataEntry(structurePath);

          std::vector<int> id = structureData.getData<int>("id");
          int N = id.size();
          if (N != pd->getNumParticles()){/*Emit error*/}

          std::vector<std::string> type = 
          structureData.getData<std::string>("type");

          std::vector<int> resBuffer;
          if(structureData.isDataAdded("resId")){
              resBuffer = structureData.getData<int>("resId");
          } else {
              resBuffer.resize(N);
              std::fill(resBuffer.begin(),resBuffer.end(),0);
          }
          // Same for chain, model and batch ...
          try {
            auto typeParamHandler = gd->getTypes();
            const int * sortedIndex = 
            pd->getIdOrderedIndices(access::location::cpu);

            auto pos  = pd->getPos(access::location::cpu,  access::mode::write);
            auto res  = pd->getResId(access::location::cpu,access::mode::write);
            ...
            for(int i=0;i<N;i++){
              int id_ = id[i];
              pos[sortedIndex[id_]].w = int(typeParamHandler->getTypeId(type[i]));
              res[sortedIndex[id_]]   = resBuffer[i];
              ...
            }
          } catch(...) {/*Emit error*/}
          {
              //Load types info (mass,radius,...) into ParticleData
              gd->getTypes()->loadTypesIntoParticleData(pd);
          }
        }
        ...
      public:
        Topology(...):...{...}
        ...
    };

At the end of the function, the method ``loadTypesIntoParticleData`` is
called. This method updates the data for each particle based on the
types that have been defined. Once this process is complete, the
initialization of ‘ParticleData‘ is considered finished, as all data
meant to be added have been introduced either in ‘State‘ or here, in
‘Structure‘.

Force Field
^^^^^^^^^^^

The ‘Force Field‘ section in the input file serves not only as a
grouping for the simulation’s interactions but also embodies a class:

.. code-block:: cpp

    class ForceField : public Interactor {
      private:
        struct schedule{
            ullint start;
            ullint end;
            bool state; //true: on, false: off
        };
        ...
        std::map<std::string, std::shared_ptr<Interactor>> interactors;
        std::map<std::string, std::shared_ptr<Interactor>> idleInteractors;
        std::map<std::string,schedule> scheduledInteractors;

        void stopInteractor(std::string interactorName){
            // Move interactor from interactors to idleInteractors
        }
        void resumeInteractor(std::string interactorName){
            // Move interactor from idleInteractors to interactors
        }
      public:
        ForceFieldBase(std::shared_ptr<Topology>   top,
                       std::string name):Interactor(top->getParticleGroup(),name),
                                         top(top){}
            ...
            interactors = top->getInteractors();
            ...
            for(auto &interactor: interactors){
                //Read interactor parameters
                DataEntry data = ...
                if(data.isParameterAdded("endStep") or
                   data.isParameterAdded("startStep")){
                    schedule sched;
                    sched.start = data.getParameter<ullint>("startStep",0);
                    ...
                    scheduledInteractors[interactor.first] = sched;
                }
            }
            ullint step = this->gd->getFundamental()->getCurrentStep();
            // Init scheduledInteractors state according current step
            ...
        }
        void sum(Computables comp,cudaStream_t st) {
            ullint step = this->gd->getFundamental()->getCurrentStep();
            for(auto &interactor: scheduledInteractors){
                //Check if interactor is active and move to/from idleInteractors
            }
            for(auto &interactor: interactors){
                interactor.second->sum(comp,st);
            }
        }
    };

The concept behind ‘ForceField‘ is to create an interactor that
encompasses all interactions. This arrangement facilitates the execution
of common operations, simplifications, or the addition of shared
features. Currently, ‘ForceField‘ is primarily used for the latter
purpose. ‘ForceField‘ analyzes the added interactors, checks if they
have defined parameters like ‘startStep‘ or ‘endStep‘, and if so, adds
them to the ‘scheduledInteractors‘ dictionary.

Upon declaration, ‘ForceField‘ requests all interactors from ‘Topology‘
(those declared within the ‘Force Field‘ section), storing them in the
‘interactors‘ dictionary, which maps each interactor (the object
instance) to a string, its name. When evaluating the ‘ForceField‘
interactors, it checks which interactors are active. These include those
not considered scheduled or, among the scheduled, those for which
‘currentStep‘ is greater than ‘startStep‘ and less than ‘endStep‘, hence
active. Active interactors are ensured to be in the list of interactors
to be computed, while inactive ones are placed in the ‘idleInteractors‘
list.



.. [1]
   As observed, we do not use UAMMD’s own ``ParticleData``, but rather
   ``ExtendedParticleData``, a class derived from ``ParticleData`` used
   by UAMMD-structured. The main difference is that this class processes
   the input to load the data present in *State*

.. [2]
   Some details have been omitted. For example, everything related to
   restarts. This process is a particular case and will be discussed in
   the future.
