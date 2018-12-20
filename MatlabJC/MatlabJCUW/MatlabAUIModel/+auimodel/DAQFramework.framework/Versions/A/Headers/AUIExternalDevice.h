//
//  AUIExternalDevice.h
//  AcqUI
//
//  Created by Barry Wark on 9/28/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


/*
 * Protocol for objects that describe external devices. Examples would include
 * amplifiers (one external device per channel), LEDs, video monitors, etc.
 *
 * Devices are characterized by having a single each for input and output and 
 * a singlemode.
 *
 * Gain is used by the AUIDAQStreamInfo to scale the output/input appropriately.
 * Output gain refers to signals leaving DAQ interface (input to the device).
 * Input gain refers to the signals arriving at the DAQ interface (output from the device).
 *
 * An implementing object will have to be coded for each device on the rig.
 *
 * All conforming objects must also implement NSCoding.
 */

@protocol AUIExternalDevice <NSObject, NSCoding>
+ (NSString*)deviceTypeName;
- (NSString*)name;
- (NSView*)configView; //if nil, no config view will be shown
- (double)outputGain; //DAC
- (double)inputGain; //ADC
- (double)mode;

// SI units for input values (e.g., "pA" for pico amps). Data in the database will be stored as values with these units.
- (NSString*)unitsForInput;
- (NSString*)unitsForInputWithMode:(int)mode;

// SI units for output values (e.g., "V" for volts). The output gain will convert Volts to this unit for output through the device.
- (NSString*)unitsForOutput;
- (NSString*)unitsForOutputWithMode:(int)mode;

- (NSString*)uuid; //UUID for this device instance
+ (id)device; //returns a new (autoreleased) instance.
@end
