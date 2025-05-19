%% FiniteQueue class
classdef FiniteQueue < handle
    properties (SetAccess = 'private', GetAccess = 'private')
        storage;
        capacity;
        ini;
        fin;
    end
    methods

        %% constructor
        function obj = FiniteQueue(capacity_,storage_)
            if nargin == 1 && capacity_ > 0
                obj.capacity = capacity_;
                obj.storage = zeros(1,capacity_);
                obj.reset();
            elseif nargin == 2 && capacity_ > 0 
                storage_size = size(storage_);
                if storage_size(1) == 1
                    obj.capacity = capacity_;
                    obj.ini = 1;
                    if storage_size(2) >= capacity_
                        obj.fin = capacity_;
                        obj.storage = storage_(1:capacity_);
                    else
                        obj.fin = length(storage);
                        obj.storage = [storage_ ...
                            zeros(1,capacity_-obj.fin)];
                    end
                end
            else
                error('Invalid input arguments.')
            end
        end

        %% makes the queue empty
        function [] = reset(obj)
            obj.ini = 1;
            obj.fin = 0;
        end

        %% destructor
        function [] = del(obj)
            obj.reset();
            delete(obj);
        end

        %% returns the queue's capacity
        function value = getCapacity(obj)
            value = obj.capacity;
        end

        %% tells the current length of the queue
        function L = length(obj)
            if obj.isEmpty()
                L = 0;
                return
            end
            L = obj.fin-obj.ini+1;
            if obj.fin < obj.ini
                L = obj.capacity-obj.ini+1+obj.fin;
            end
        end

        %% returns the full queue as a row vector
        function [queue,L] = q2vec(obj)
            L = obj.length();
            if L == 0
                queue = [];
                return
            end
            queue = zeros(1,L);
            i = obj.ini-1;
            for i_ = 1:L
                i = i + 1;
                if i > obj.capacity
                    i = 1;
                end
                queue(i_) = obj.storage(i);
            end
        end

        %% doubles the queue size
        function [] = resize(obj)
            [queue,L] = obj.q2vec();
            obj.capacity = obj.capacity * 2;
            obj.storage = [queue zeros(1,obj.capacity-L)];
            obj.ini = 1;
            obj.fin = L;
        end

        %% tells whether the queue is empty
        function bool_val = isEmpty(obj)
            bool_val = (obj.fin == 0 && obj.ini == 1);
        end

        %% enqueues a new element new_el
        function [] = enqueue(obj,new_el)
            if obj.length() == obj.capacity
                obj.resize();
            end
            if obj.fin >= obj.ini && obj.fin < obj.capacity || ...
                    obj.fin < obj.ini-1 || obj.isEmpty
                obj.fin = obj.fin + 1;
                obj.storage(obj.fin) = new_el;
            elseif obj.ini > 1 && obj.fin == obj.capacity
                obj.fin = 1;
                obj.storage(obj.fin) = new_el;
            end
        end

        %% dequeues the first element enqueued
        function [] = dequeue(obj)
            if obj.isEmpty
                return
            end
            if obj.ini == obj.fin
                obj.reset();
                return
            end
            obj.ini = obj.ini + 1;
            if obj.ini > obj.capacity
                obj.ini = 1;
            end
        end

        %% returns the first element enqueued
        function el = front(obj)
            if obj.isEmpty
                el = [];
                return
            end
            el = obj.storage(obj.ini);
        end

    end
end