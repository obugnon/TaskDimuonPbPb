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
                                TaskDimuonPbPb(const char *name);
        virtual                 ~TaskDimuonPbPb();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list
        TH1F*                   fHistPt;        //! dummy histogram

        TaskDimuonPbPb(const TaskDimuonPbPb&); // not implemented
        TaskDimuonPbPb& operator=(const TaskDimuonPbPb&); // not implemented

        ClassDef(TaskDimuonPbPb, 1);
};

#endif
