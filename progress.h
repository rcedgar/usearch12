#pragma once

class SeqSource;

void StopProgressThread();
void StartProgressThread();

typedef void (*PTR_PROGRESS_CB)(string &s);

void ProgressNote(const char *fmt, ...);
void ProgressNoteNoPrefix(const char *fmt, ...);
void ProgressNoteLog(const char *fmt, ...);

void ProgressStartOther(const string &Str,
  PTR_PROGRESS_CB CB = 0);
void ProgressDoneOther();

void ProgressStartSS(SeqSource &SS, const string &Str,
  PTR_PROGRESS_CB CB = 0);
void ProgressDoneSS();

uint32 *ProgressStartLoop(uint32 N, const string &Str,
  PTR_PROGRESS_CB CB = 0);
void ProgressDoneLoop();

void SearcherCB(string &str);
void ClusterCB(string &str);
void ClosedRefCB(string &str);
void OtuTabCB(string &str);
void UPARSECB(string &str);
void SintaxCB(string &str);

#define basenm(fn)	string(BaseName(fn.c_str()))
