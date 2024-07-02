#pragma once

#include"InputOutput/Input/Input.cuh"

//Units
#include"GlobalData/Units/UnitsHandler.cuh"
#include"GlobalData/Units/UnitsLoaders.cuh"

//Fundamental
#include"GlobalData/Fundamental/FundamentalHandler.cuh"
#include"GlobalData/Fundamental/FundamentalLoaders.cuh"

//Ensemble
#include"GlobalData/Ensemble/EnsembleHandler.cuh"
#include"GlobalData/Ensemble/EnsembleLoaders.cuh"

//Types
#include"GlobalData/Types/TypesHandler.cuh"
#include"GlobalData/Types/TypesLoaders.cuh"


namespace uammd{
namespace structured{

    class GlobalDataBase{

        private:

            std::shared_ptr<ExtendedSystem> sys;

            //////////////////////////////////////////

            std::shared_ptr<Units::UnitsHandler>       unitsHandler;
            std::shared_ptr<Ensemble::EnsembleHandler> ensembleHandler;

            std::vector<std::string> path;

            //////////////////////////////////////////

        public:


            GlobalDataBase(std::shared_ptr<ExtendedSystem> sys,
                           std::vector<std::string> path);

            std::shared_ptr<ExtendedSystem> getSystem(){return sys;}

            std::shared_ptr<Units::UnitsHandler>             getUnits(){return unitsHandler;}
            std::shared_ptr<Ensemble::EnsembleHandler>       getEnsemble(){return ensembleHandler;}

            void updateInputGlobalData();
    };

    class GlobalData{

       private:

            std::shared_ptr<GlobalDataBase>    globalDataBase;
            std::vector<std::string>           path;

            std::shared_ptr<Types::TypesHandler>             typesHandler;
            std::shared_ptr<Fundamental::FundamentalHandler> fundamentalHandler;

            void init();

        public:

            GlobalData(std::shared_ptr<GlobalDataBase>  globalDataBase,
                       std::vector<std::string>         path);

            GlobalData(std::shared_ptr<ExtendedSystem>  sys,
                       std::vector<std::string>         path);

            GlobalData(std::shared_ptr<ExtendedSystem>  sys);

            std::shared_ptr<GlobalDataBase> getGlobalDataBase(){return globalDataBase;}

            std::shared_ptr<ExtendedSystem> getSystem(){ return this->getGlobalDataBase()->getSystem();}

            std::shared_ptr<Units::UnitsHandler>             getUnits()      {return this->getGlobalDataBase()->getUnits();}
            std::shared_ptr<Ensemble::EnsembleHandler>       getEnsemble()   {return this->getGlobalDataBase()->getEnsemble();}
            std::shared_ptr<Types::TypesHandler>             getTypes()      {return typesHandler;}
            std::shared_ptr<Fundamental::FundamentalHandler> getFundamental(){return fundamentalHandler;}

            void updateInputGlobalData();
    };

}}
