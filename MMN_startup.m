   
% start up code for active x control of the RP2
    clear;
    close all;
    % set up actx server/control 
    handles.RP = actxcontrol('RPco.x');
    RP	= handles.RP;

    % connect to the device and halt any ongoing processes
    RP.ConnectRP2('USB',1);
    RP.Halt;
    RP.ClearCOF;
    
    % load our rcx file and run it
    %RP.LoadCOF('C:\TDT\RPvdsEx\Examples\TDT_test_soft_trig.rcx');
    RP.LoadCOF('C:\Users\seema\OneDrive\Documents\MATLAB\MMN.rcx');
    RP.Run; %
    %close all;