
function [assigned_label, right_or_wrong] = KNN(training_data , training_labels, k, distance, test_object, test_label, graphical)
% Assigned label: el label que t'assigna.
% right_or_wrong: et diu si esta be comparantho amb el label que li passes.

%training_data:each column represents a training object; and the rows of
% each column represent the attributes corresponding to a training object.

%training_labels: this parameter is a row vector with N elements, all of them integer values. 
%Each element in training_labels works as an identification label of the category of the corresponding
%training object in the training_data matrix. That is, the 1st element of training_labels is 
%the category label of the training object contained in the 1st column of the training_data matrix, and so on.

%k: this parameter is an integer indicating the number of nearest neighbours that have to be employed by the classification algorithm.
%Its minimum value is 1.

%distance: this parameter is a text string that indicates the type of distance employed by the k-NN classifier. 
%It can take one of the following values:
%?euclidean?: Euclidean distance
%?manhattan?: Manhattan distance
%?chebysev?: Chebysev distance


%test_object: this parameter is a column vector with M elements, and it
%represents the object we want to classifiy. Notice that it is a single vector, so it
%represents a single object.

%test_label: this parameter is an integer that indicates the category label that
%corresponds to the test object. This parameter must only be employed for
%evaluating if the classifier made a right or wrong classification decision.

%graphical: this parameter is binary. Its goal is to provide a visual representation of the classification process.
%If it equals 1, the function must display the objects in the training data as points in space, using a different color 
%for each category (which is specified by the training_labels vector).
%Moreover, it also must display the test object by a point of a different color. If this parameter equals 0, the
%function displays nothing.

    training_objs = size(training_data,2);
    obj_comp = size(training_data,1);

    D = zeros(1,training_objs);
    %Nem a mirar per a totes les dades de training. 
    for i=1:training_objs
        val = 0;
        if(strcmp('manhattan',distance))
            %Haurem de fer la simple resta entre els atributs.
            for j=1:obj_comp
            	aux = abs(training_data(j,i)-test_object(j));
                val = val + aux;
            end
            
            D(i) = val;
        end
        if(strcmp('euclidean',distance))
            
            for j=1:obj_comp
            	aux = (training_data(j,i)-test_object(j))^2;
                val = val + aux;
            end
            
            D(i) = sqrt(val);
        end
        if(strcmp('chebyshev',distance))
            
            for j=1:obj_comp
            	val(j) = abs(training_data(j,i)-test_object(j));
            end
            
            [D(i),pos] = max(val);
        end
        
    end
    %ara que ja tenim fetes les operacions (distancies) mirem quines son les k minimes.
    
    for n=1:k
        [value,pos] = min(D);
        %mirem quina label es.
        Result(n) = training_labels(pos);
        %després posem la D a infinit perque no surti repetida.
        D(pos) = inf;
    end
    %el label guanyador sera aquell que s'hagui repetit més vegades.
    [unique_strings, ~, string_map] = unique(Result);
    assigned_label = unique_strings(mode(string_map));
    right_or_wrong = 0;
    
    if ( strcmp(assigned_label,test_label) )
        right_or_wrong = 1;
    end
    
    if (graphical == 1)
        figure,  
        scatter(test_object(1), test_object(2),150,'+','m','LineWidth',2);
        hold on;
        gscatter(training_data(1,:), training_data(2,:),training_labels,'rgbm', 'xo*-');
    end
end
