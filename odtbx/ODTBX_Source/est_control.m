classdef est_control < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Estimator
        est_type
        myest
        
        
    end
    
    methods
        
        function obj = est_control(varargin)
            %% Input Parsing and Setup            
            if nargin >= 1,
                obj.est_type = varargin{1};
                obj.myest = feval(obj.est_type, varargin{2:end});
            else
                disp "Must provide inputs!";
            end
        end
        
        
        function set_controllers(obj,varargin)
            % This function needs to be run to substitute the
            % events/controls functions from this class (or other specified functions)
            % for the defaults located in the estimator.
            if nargin >= 2,
                obj.myest.events_fcn = varargin{1};
                obj.myest.control_events_fcn = varargin{2};
            else
                obj.myest.events_fcn = @obj.events;
                obj.myest.control_events_fcn = @obj.control_events;
            end
        end
        
        
        function varargout = run_sim(obj)
            %% Run estimator
            [t,Xhat,Phat,e,Y] = obj.myest.run_estimator();
            
            %% Output results
            if nargout >= 3,
                varargout{1} = t;
                varargout{2} = Xhat;
                varargout{3} = Phat;
            end
            if nargout >= 4,
                varargout{4} = e;
            end
            if nargout >= 5,
                varargout{5} = Y;
            end
        end
        
        
        function [value,isterminal,direction] = events(obj,t,X,varargin)
            % This function is used to kick the integrator out of its loop
            % at a certain point (as determined by time or state
            % conditions) in order to perform an action desginated by
            % control_events.
            
            % See header in integev.m for details on event function formats
            
%             disp "Controller"
            % Consider using functions for conditions
            
            %% Event 1:
            condition1 = X(1); %  We can change this to be anything related to t or X
            terminal1 = 0;
            direction1 = 0;
            
            %% Event 2:
            condition2 = t - 295; %  We can change this to be anything related to t or X
            terminal2 = 1;
            direction2 = 0;
                        
            %% Put all the values together to be returned
            value = [condition1; condition2]; % Condition to look for
            isterminal = [terminal1; terminal2]; % Do we need to end the integration at this point?
            direction = [direction1; direction2]; % Is there a direction involved?
        end % events
        
        
        function [X_state_mod, P_mod] = control_events_default(obj,t,X,P,varargin)
            
            % This function is used to change the state/covariance once a
            % condition has been detected.
            X_state_mod = X;
            P_mod = P;
            
        end % control_events
        
    end % methods
    
end % est_control

