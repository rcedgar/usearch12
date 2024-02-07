#pragma once

class SeqSource;

void StopProgressThread();
void StartProgressThread();

typedef void (*PTR_PROGRESS_CB)(string &s);

void ProgressNote(const char *fmt, ...);
void ProgressNoteLog(const char *fmt, ...);

void ProgressStartOther(const char *fmt, ...);
void ProgressDoneOther();

void ProgressStartCB(PTR_PROGRESS_CB CB, const char *fmt, ...);
void ProgressDoneCB();

void ProgressStartSS(SeqSource &SS, const char *fmt, ...);
void ProgressDoneSS();

uint32 *ProgressStartLoop(uint32 N, const char *fmt, ...);
void ProgressDoneLoop();
