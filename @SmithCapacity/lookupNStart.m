function nStart = lookupNStart(A)
% lookupNStart returns an appropriate nStart for a given A (based on
% previous trials).

% table columns are [A, nStart]
table = [ ...
    32, 43;
    28, 38; % Up to A=30
    27, 34;
    26, 33;
    25, 32;
    24, 30;
    23, 29;
    22, 28;
    21, 27;
    20, 25;
    19.3, 27;
    18.5, 26;
    17.9, 25;
    17.7, 23;
    16.9, 22;
    16.1, 21;
    15.3, 20;
    14.5, 19;
    13.7, 18;
    12.9, 17;
    12.1, 16;
    11.3, 13;
    10.5, 12;
    9.7, 11;
    8.8, 10;
    7.9, 9;
    7.0, 8;
    6.1, 7;
    5.1, 6;
    4.1, 5;
    3.0, 4;
    1.7, 3;
    0, 2];

indx = find(table(:,1) <= A, 1);

nStart = table(indx, 2);

end

