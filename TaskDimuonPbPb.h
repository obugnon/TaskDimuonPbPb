/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef TaskDimuonPbPb_H
#define TaskDimuonPbPb_H

#include "AliAnalysisTaskSE.h"

class TaskDimuonPbPb : public AliAnalysisTaskSE  
{
    public:
                                TaskDimuonPbPb();
                                TaskDimuonPbPb(const char *name,int firstRun, int lastRun, UInt_t triggerClass);
        virtual                 ~TaskDimuonPbPb();
        
        virtual void            NotifyRun();
        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;           //! input event
        AliMuonTrackCuts*       fMuonTrackCuts; //! usual cuts on single muon tracks
        UInt_t                  fTriggerClass;  //! trigger selection
        int                     fFirstRun, fLastRun;

        TH1I *fHistoTotalEventsPerRun;      //! histogram to store number of events
        TH1I *fHistoPSEventsPerRun;         //! for the physics selection
        TH1I *fHistoEventsBeforePSPerRun;   
      
        TH1I *fHistoCMULEventsInCINT7;      //! for the normalization
        TH1I *fHistoCMSLEventsInCINT7;
        TH1I *fHistoCMULEventsInCMSL;
      
        TH1F *fHistoNumberMuonsCuts;

        THnSparseD *fHistoDiMuonOS;     //! histogram to store some properties of dimuons unlike signe
        //THnSparseD *fHistoDiMuonLS;
      
        THnSparseD *fHistoSingleMuon;   //! histogram to store some properties of single muons
      
        TList *fListEventHistos;        //! list to save the events histograms 
        TList *fListSingleMuonHistos;   //! list to save the single muon histograms
        TList *fListDiMuonHistos;       //! list to save the dimuon histograms

        TaskDimuonPbPb(const TaskDimuonPbPb&); // not implemented
        TaskDimuonPbPb& operator=(const TaskDimuonPbPb&); // not implemented

        ClassDef(TaskDimuonPbPb, 2);
};

#endif
