%design for D(x)
function D = stiffnessD(question_def, E, v)
%strain
if question_def == 1
    D=E/(v+1)/(2*v-1) *[v-1, -v, 0; -v, v-1, 0; 0, 0, 0.5*(2*v-1)];
else
    %stress
    if question_def == 2
        D=1/E *[1, -v, 0; -v, 1, 0; 0, 0, 2*(1+v)];
    else
        error("only for 1 and 2")
    end
end
end


