/*
 *  UserTrigger.h
 *  Acquirino
 *
 *  Created by CJ Bell on 12/5/06.
 *
 *  Returns true when the user has "pressed the trigger"
 *
 *  A "trigger" is represented by a message-box with a single button.
 */
 
#ifndef INCLUDED_UserTrigger_h
#define INCLUDED_UserTrigger_h

/**
 * If a "trigger" message-box is not already presented to the user, then this
 * will create one. When the user clears the message-box, the next call to this
 * function will return true. Otherwise, this will return false.
*/
int Triggered();

#endif // INCLUDED_UserTrigger_h
