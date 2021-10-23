#pragma once


typedef uint16_t node_t; 

const node_t cubic_neighbours[512*3]  = {76,180,72,86,388,79,89,65,271,95,66,91,148,267,272,81,262,480,113,114,111,127,83,121,92,318,137,138,85,256,108,97,244,155,482,181,103,124,158,160,100,336,431,104,436,505,449,511,495,112,503,499,510,486,122,116,334,171,123,489,487,163,235,131,129,284,485,221,366,143,134,415,465,232,471,448,383,455,290,164,460,412,509,322,347,424,369,493,255,421,329,490,227,405,393,410,475,427,402,468,306,234,464,467,172,389,472,274,459,492,381,166,210,385,74,419,107,313,315,170,288,314,242,443,350,461,430,434,191,426,429,195,291,442,488,473,253,201,394,280,354,345,399,340,351,441,276,338,167,222,179,231,379,265,173,408,416,463,428,186,182,323,297,250,362,339,197,398,356,189,275,299,216,230,263,332,346,380,266,409,213,303,273,316,376,326,287,270,198,259,304,285,414,149,327,262,2,66,65,3,159,392,420,226,70,451,69,68,165,271,81,484,68,388,451,469,94,0,74,333,470,352,72,38,114,114,445,94,79,0,135,237,375,501,343,169,378,479,1,76,184,305,178,271,5,70,83,208,479,135,7,82,318,496,367,203,9,318,204,1,87,86,496,92,92,413,101,91,2,90,89,165,246,93,3,89,87,8,88,102,142,91,75,192,72,239,3,115,109,118,106,482,10,239,219,205,477,272,282,108,108,13,124,88,165,204,246,126,93,148,12,278,449,14,258,480,390,431,217,96,140,111,38,109,99,10,100,107,96,321,446,452,456,375,6,107,321,16,375,188,6,237,74,6,75,95,142,506,117,18,188,510,400,116,469,96,419,121,192,196,156,324,365,123,7,119,196,18,163,163,19,121,100,12,272,301,208,127,144,102,141,125,7,171,224,422,348,221,21,301,256,360,131,130,21,132,131,476,136,146,225,268,136,23,256,76,192,83,132,139,134,218,8,138,137,9,143,403,136,478,106,484,480,225,126,218,115,93,466,138,23,225,466,126,146,149,200,466,144,133,232,207,474,438,258,4,103,232,64,145,181,200,414,498,261,418,414,211,154,373,264,175,152,206,324,244,11,156,155,120,383,433,504,435,162,12,160,66,282,267,158,13,174,483,319,397,164,277,158,122,20,123,292,26,162,90,69,101,381,37,417,417,49,264,220,183,500,328,78,303,190,39,288,127,19,221,385,34,234,373,51,231,160,286,292,408,153,338,418,228,426,349,382,444,182,80,456,363,50,250,419,0,388,324,11,150,250,53,178,384,508,168,461,80,186,233,371,461,184,53,197,309,202,302,116,445,113,289,56,194,231,283,170,195,42,291,135,94,119,308,359,457,189,453,488,348,43,191,119,445,122,186,55,233,259,62,302,492,211,327,150,506,145,260,45,209,330,187,352,367,360,85,101,451,86,394,98,474,458,154,342,392,147,391,367,82,125,201,491,224,464,37,492,462,152,199,296,320,368,265,60,377,500,391,243,473,293,332,332,57,340,431,295,106,141,413,137,438,98,251,222,386,168,171,22,129,358,49,220,251,401,266,209,128,296,143,133,141,474,67,397,457,30,405,176,387,432,383,286,336,439,57,263,173,50,190,146,24,149,197,407,185,172,33,389,490,20,334,255,400,499,113,77,510,263,502,454,97,282,95,378,374,409,427,476,284,264,40,381,266,214,438,336,10,155,460,252,443,90,413,102,312,254,273,340,437,473,297,254,361,179,54,182,219,335,223,245,292,447,364,45,406,249,377,247,463,29,236,134,9,130,300,396,376,104,390,148,310,63,198,470,293,201,387,151,486,372,5,65,230,58,238,167,153,242,379,51,213,223,59,243,159,4,372,471,133,415,294,502,287,309,62,279,69,2,81,124,4,99,247,60,343,508,35,402,311,56,341,488,48,411,337,278,162,103,277,449,270,502,346,399,46,457,361,298,398,99,159,239,363,344,190,241,21,485,287,63,326,174,229,450,269,62,285,170,40,373,411,407,189,337,26,382,191,44,507,174,252,164,215,260,307,326,396,269,321,217,395,224,212,470,323,54,249,370,281,331,345,57,317,409,401,257,129,360,125,198,187,357,169,60,408,316,63,312,446,80,350,472,33,478,346,293,333,397,193,394,333,187,270,331,353,259,370,355,275,304,353,247,458,39,344,459,40,315,314,39,342,343,61,304,299,335,477,85,8,84,359,161,366,507,212,348,109,295,112,369,27,433,398,53,297,154,120,181,368,453,356,285,61,294,199,64,465,338,374,169,487,30,359,356,355,202,357,298,310,215,58,216,307,73,309,235,18,494,251,317,439,229,13,244,511,277,290,175,49,328,341,55,370,216,47,248,275,407,339,315,206,462,273,78,316,313,283,452,354,47,299,279,58,307,495,28,387,320,128,195,435,177,351,305,41,447,349,48,442,202,73,368,310,361,312,477,46,345,330,311,357,325,56,330,302,355,331,380,374,222,329,193,319,130,203,301,249,353,281,377,54,379,456,283,179,416,491,253,455,120,458,319,22,489,84,208,203,352,212,325,432,28,322,339,298,311,411,440,185,267,390,262,288,153,173,328,358,240,112,77,111,257,61,378,213,254,362,376,78,240,362,50,265,500,59,358,242,36,166,290,177,504,156,25,229,389,183,497,497,37,172,417,497,220,347,228,261,180,1,71,234,35,384,105,372,258,207,214,508,402,67,207,493,31,490,308,46,205,503,295,436,257,454,294,226,161,308,281,55,323,404,47,280,494,117,236,300,223,423,274,32,392,415,139,468,405,437,399,227,31,404,253,437,410,289,233,341,303,51,175,240,59,300,406,31,421,276,371,289,505,27,424,246,88,218,150,64,152,268,23,403,421,52,364,166,386,167,428,151,176,118,38,180,483,67,427,410,29,416,426,128,425,454,401,439,412,28,503,422,491,428,176,43,422,420,32,241,425,52,418,430,43,432,433,42,429,105,14,217,429,228,369,322,157,430,442,42,435,434,157,349,395,14,505,404,406,248,243,147,219,423,335,230,443,371,441,440,48,444,351,44,434,245,41,440,441,177,460,196,75,188,448,110,305,350,252,450,450,25,446,278,15,104,447,286,448,204,68,71,344,110,455,507,194,325,238,396,423,452,25,365,178,110,363,280,193,227,365,206,313,462,36,314,444,26,245,185,41,184,342,211,459,498,52,255,465,34,210,327,24,464,145,142,144,468,34,471,403,33,467,71,484,118,296,73,260,467,24,268,475,35,306,248,45,215,205,147,226,481,32,472,481,132,241,317,98,354,306,139,481,82,496,79,140,5,105,478,476,475,506,11,97,485,161,420,469,70,140,284,22,483,261,17,501,489,20,329,194,44,276,366,19,487,393,30,235,364,425,209,210,36,199,494,29,393,334,400,493,501,16,347,87,479,84,384,386,385,499,151,463,236,17,498,168,214,380,486,77,495,238,279,269,424,16,395,382,157,509,436,15,412,115,200,482,291,453,320,391,183,274,504,27,511,237,17,117,509,15,337};

const node_t next_on_face[512*3] = {79,419,94,204,180,479,91,262,69,239,65,93,258,159,124,271,372,140,188,74,375,125,135,123,87,85,218,137,203,134,99,482,336,244,506,324,148,100,162,158,108,229,105,449,395,436,278,509,501,321,424,236,237,261,196,117,235,127,163,366,489,122,490,130,221,241,284,171,319,138,136,268,327,146,467,450,156,452,337,292,444,505,504,369,495,412,432,494,463,410,487,393,457,227,493,406,481,420,274,403,472,172,465,468,385,234,475,508,462,210,242,381,464,497,72,118,111,458,314,190,170,459,264,245,305,185,433,442,195,176,430,348,191,351,194,248,364,260,308,399,477,354,404,216,349,440,488,175,417,358,363,173,362,379,373,303,421,498,425,184,250,398,323,179,377,341,186,281,325,289,311,345,332,439,230,215,279,500,223,240,265,169,247,343,257,285,269,309,259,310,316,287,150,232,199,5,89,159,2,95,267,402,483,474,484,204,271,451,90,81,5,469,69,1,68,118,192,180,114,307,296,202,0,107,75,6,196,72,1,72,83,113,112,486,273,328,376,496,388,135,461,446,182,2,480,68,7,367,79,192,127,479,8,479,203,360,138,84,451,79,92,1,84,88,8,246,204,3,271,246,2,101,102,142,66,90,496,137,101,126,115,89,445,135,74,282,91,506,107,469,217,11,108,95,438,394,317,4,239,100,10,160,272,413,69,86,413,144,91,4,158,449,15,431,148,5,258,217,295,118,480,6,419,321,282,244,124,38,106,112,448,344,178,77,114,109,295,495,111,445,111,510,38,113,94,3,466,482,400,122,113,17,494,188,484,109,180,7,94,122,155,154,455,19,83,196,445,334,123,20,171,119,13,103,99,360,82,171,466,246,225,208,121,221,209,426,320,22,131,125,9,301,132,360,284,136,21,481,134,144,143,471,139,143,130,0,119,82,476,403,256,413,318,143,8,256,225,415,132,306,96,70,105,133,102,137,95,102,145,9,415,141,142,141,232,64,506,144,126,268,149,392,205,243,390,272,278,24,414,466,11,145,152,499,387,428,64,462,324,288,167,408,211,458,181,10,181,383,11,365,229,322,382,434,277,124,174,3,99,372,12,336,292,485,359,226,26,278,160,18,487,121,252,290,158,89,68,88,36,385,167,386,338,242,386,508,380,374,343,408,283,315,373,7,489,129,37,467,389,153,265,190,13,450,164,51,264,328,151,432,422,435,290,441,53,305,363,283,379,182,38,76,71,120,482,414,54,186,456,389,391,220,41,178,197,407,440,184,80,323,233,333,330,198,18,75,237,407,356,488,50,344,288,43,434,507,76,75,121,397,329,280,56,507,276,128,429,291,192,188,163,53,339,185,63,270,357,36,152,465,181,115,149,293,253,224,355,309,368,208,130,318,165,71,87,46,219,226,365,152,315,67,438,508,84,83,301,45,425,296,34,166,199,342,414,492,224,507,352,51,273,362,168,207,266,45,307,216,58,299,248,14,321,140,126,88,138,147,477,223,49,497,500,19,485,301,374,167,168,335,300,243,491,348,470,23,146,218,147,420,308,193,490,404,418,347,429,25,174,244,335,216,238,51,179,170,133,465,145,55,289,461,34,306,384,30,163,494,29,117,498,6,501,117,58,269,423,10,159,115,78,358,300,32,132,485,153,314,166,59,391,219,13,97,156,26,447,440,165,218,93,353,377,343,47,406,215,54,247,281,50,297,178,98,439,266,460,174,350,491,473,410,297,213,312,52,493,499,23,85,131,401,294,378,14,372,103,353,285,302,73,215,209,228,498,501,390,81,66,57,346,454,49,373,381,50,408,377,401,380,438,282,148,262,24,225,403,396,279,285,187,287,346,165,65,70,12,267,108,254,303,316,183,472,392,355,189,339,44,441,289,511,103,164,12,337,104,62,238,307,47,394,227,353,370,323,272,66,97,456,313,231,476,129,483,62,304,294,160,383,447,502,198,326,39,242,173,371,341,194,277,460,504,42,488,320,286,245,162,473,470,346,61,454,287,109,431,503,128,368,260,53,362,361,339,361,357,47,230,477,59,423,376,21,203,127,62,202,331,78,213,175,61,259,247,110,184,447,35,468,481,58,260,309,161,457,205,73,302,279,298,312,198,298,330,341,63,361,273,206,170,452,36,288,342,40,313,462,78,326,312,57,251,354,9,92,367,193,483,489,453,296,195,96,395,375,28,509,430,55,182,249,206,156,150,212,194,330,63,376,269,211,149,464,49,240,303,20,227,319,56,357,352,355,281,259,293,263,340,293,352,270,20,116,493,219,299,423,286,100,155,15,162,382,153,222,169,407,398,311,57,399,473,56,233,370,39,154,459,60,378,304,39,363,455,46,340,317,502,332,333,16,369,261,212,422,191,157,444,442,80,443,450,177,276,434,187,470,325,331,249,304,98,280,299,356,370,302,453,275,202,187,311,310,59,328,220,30,308,366,256,367,129,254,310,398,254,250,265,110,190,250,52,209,406,25,324,313,161,221,487,496,125,85,73,320,356,228,424,433,55,331,275,276,443,233,4,105,65,40,175,231,338,380,378,16,237,107,396,316,240,60,249,379,61,169,409,54,231,213,214,409,222,40,492,417,26,349,509,120,448,336,35,168,385,386,210,234,166,384,222,28,176,486,0,86,469,33,274,497,480,267,104,147,500,274,32,226,391,29,405,235,193,354,474,16,217,505,300,238,326,67,319,394,298,197,297,437,345,457,334,510,255,409,251,454,35,427,207,23,478,467,31,248,280,30,410,399,45,404,421,411,197,275,60,173,338,374,266,257,437,393,416,48,185,189,15,322,503,90,92,141,200,327,154,133,134,468,29,428,253,37,220,264,52,261,426,96,74,388,161,392,241,31,255,364,43,224,428,396,223,230,27,347,395,128,364,418,228,195,425,67,475,284,491,463,176,42,426,369,157,191,432,390,436,106,43,387,322,27,435,429,44,430,349,42,504,351,295,104,412,405,253,340,214,474,251,401,317,263,41,411,444,371,351,460,48,291,435,252,461,441,48,382,245,119,114,116,25,456,350,41,292,448,286,455,305,277,505,258,252,229,446,101,70,388,283,446,365,291,189,368,502,257,439,110,383,458,80,452,179,46,359,405,120,342,344,211,381,315,177,164,443,371,350,186,206,199,314,151,416,236,24,172,492,64,471,210,200,93,146,33,464,268,139,234,471,451,140,419,212,333,201,34,232,415,32,389,478,437,201,332,98,207,397,476,402,306,478,131,427,335,205,345,33,136,475,208,87,76,484,262,431,139,241,472,200,155,239,22,397,427,71,81,106,21,366,420,151,510,495,19,235,359,453,442,411,22,123,329,31,329,334,416,422,201,37,459,327,400,421,490,18,236,393,77,503,387,86,82,318,183,417,172,17,418,255,400,486,463,183,243,358,17,375,347,263,270,294,28,112,436,177,433,511,14,511,424,142,150,97,44,325,348,214,384,402,157,412,337,77,499,116,27,449,290};

const node_t prev_on_face[512*3] = {135,388,74,87,71,76,90,66,81,115,159,89,103,372,99,70,65,105,237,75,107,171,82,119,88,84,138,143,318,130,100,239,155,156,97,150,278,272,160,174,124,244,217,258,505,412,104,337,347,375,395,498,117,501,163,188,494,221,121,487,329,123,334,132,301,485,483,129,489,225,256,403,464,149,268,446,229,365,382,162,245,424,511,433,387,503,322,393,236,416,359,235,405,404,490,421,472,241,392,467,478,389,210,471,234,384,306,402,314,199,166,417,492,172,114,180,109,344,342,288,373,315,381,440,447,184,429,435,291,422,432,191,507,434,276,215,406,209,205,457,345,299,280,248,442,444,411,328,264,220,250,190,265,213,231,175,364,255,418,197,178,297,249,182,379,370,233,323,330,194,341,317,340,263,238,216,307,358,243,300,377,408,343,304,378,294,285,279,302,198,312,326,152,145,465,372,271,3,262,91,282,207,427,397,81,71,165,70,101,2,271,140,451,180,204,484,75,76,38,309,260,368,94,419,6,74,188,192,479,180,192,510,111,495,316,303,240,82,86,0,186,350,456,69,262,484,135,125,496,76,121,208,85,87,208,367,256,8,101,388,496,204,479,8,87,218,165,93,65,165,91,69,413,102,95,2,86,318,413,246,466,3,114,119,0,97,66,142,321,419,140,506,244,282,251,474,354,124,159,10,99,336,12,92,90,451,90,141,142,258,124,277,278,436,390,140,372,14,431,109,484,375,74,96,272,97,13,111,118,295,305,455,363,112,113,38,109,503,77,116,114,77,72,111,445,239,93,200,510,334,445,237,236,18,71,106,38,123,135,445,383,181,458,163,127,192,119,116,20,122,489,7,108,158,4,129,367,7,146,93,218,301,83,19,296,425,195,171,284,360,134,203,21,256,129,476,130,241,139,232,141,415,132,415,9,79,94,7,131,478,23,141,92,9,218,85,23,468,134,481,217,469,5,143,144,413,506,91,144,137,134,133,145,102,133,232,150,142,466,225,24,391,226,219,104,267,12,146,327,200,324,506,64,463,486,176,150,199,206,173,242,338,414,342,120,336,482,120,244,324,25,430,509,349,164,103,13,65,239,4,162,100,286,420,366,308,292,337,12,196,235,19,174,460,277,246,271,204,242,210,386,166,222,153,222,384,214,338,378,60,231,313,40,125,123,22,497,464,33,288,408,50,158,229,252,303,373,49,428,387,43,351,504,460,250,184,110,456,231,54,118,72,1,154,155,200,179,323,80,497,274,500,185,305,53,197,411,41,461,182,55,270,352,357,117,196,6,411,275,453,173,363,39,348,430,44,83,72,196,394,319,227,289,325,44,320,426,42,121,75,18,184,398,407,310,287,187,210,462,64,414,482,466,470,473,491,356,302,73,84,301,9,88,68,1,308,477,147,313,324,462,402,474,214,203,479,127,260,364,128,465,385,36,459,154,327,470,348,325,379,303,254,380,508,438,248,260,58,215,230,47,105,395,96,225,246,8,243,205,335,358,417,183,127,366,21,380,338,386,219,423,59,201,422,212,138,268,126,205,392,161,280,329,31,426,261,369,156,450,13,423,299,58,373,379,283,144,471,64,186,341,371,385,468,35,393,487,18,463,494,17,188,375,17,230,279,396,482,99,3,376,328,59,420,481,21,167,288,36,223,500,147,229,108,11,444,292,41,89,88,126,304,249,60,216,404,45,323,377,353,363,362,53,438,317,401,443,164,450,416,201,437,361,362,273,498,421,400,136,138,360,409,454,61,449,105,4,331,304,62,296,307,45,347,418,17,267,480,2,439,332,502,417,175,40,362,173,60,251,409,214,66,272,390,467,146,23,326,238,62,333,198,502,68,89,5,100,148,282,312,213,78,391,389,32,370,356,407,194,351,371,290,449,158,148,162,15,309,269,58,404,354,193,249,331,55,108,267,95,179,452,170,427,131,22,269,259,61,292,336,448,294,270,63,190,314,153,276,233,56,511,164,177,195,442,453,160,447,26,332,201,333,285,257,502,112,106,436,209,320,73,398,250,254,311,398,310,354,216,335,240,223,396,221,130,208,259,309,355,328,273,51,343,285,353,448,178,41,475,234,139,279,215,73,226,359,46,307,202,62,357,361,63,339,357,56,316,310,254,365,315,283,462,242,39,459,170,206,273,376,63,345,439,98,203,137,496,329,397,22,291,368,128,107,217,16,432,412,157,281,186,54,152,365,11,352,507,56,287,316,396,492,414,24,175,358,78,489,490,193,325,311,187,302,370,353,473,346,57,346,470,187,490,122,400,223,477,230,383,160,10,509,278,26,408,167,374,275,197,298,332,345,437,311,289,55,314,458,211,247,169,61,458,190,110,477,399,57,270,263,293,501,424,228,507,224,43,434,382,48,446,461,252,435,441,44,330,333,212,259,281,247,317,394,47,202,275,331,368,189,355,198,330,298,500,240,49,487,457,161,131,85,125,297,312,298,213,297,50,178,344,50,421,425,45,452,156,206,359,485,19,318,82,360,202,296,453,429,347,27,341,281,355,289,441,461,159,258,5,170,264,51,169,222,409,321,501,6,300,326,78,265,247,54,257,343,374,377,179,51,168,266,374,264,459,37,337,444,157,155,455,286,234,508,386,384,166,34,167,385,168,495,432,151,419,79,451,172,472,183,431,262,148,392,243,183,274,420,147,494,410,30,397,280,98,424,321,14,376,423,269,474,483,193,361,339,53,405,340,46,493,116,499,257,266,439,508,475,67,268,136,33,227,406,47,457,393,437,364,248,31,189,185,339,169,265,153,378,380,401,253,405,29,488,440,407,436,509,28,102,101,137,181,149,211,471,143,139,410,463,491,381,497,49,425,498,228,469,107,0,485,226,32,406,493,52,176,348,491,238,300,335,505,369,16,426,209,52,418,429,128,483,402,476,422,416,151,433,195,228,322,434,43,480,104,295,430,176,28,369,504,42,351,191,157,442,433,177,503,431,15,399,410,473,266,207,98,454,251,57,245,185,48,443,276,177,349,488,42,460,350,371,440,349,26,122,94,113,450,452,80,305,245,286,447,383,110,103,511,14,350,174,25,86,69,469,313,456,25,320,488,356,263,294,401,344,448,120,182,446,283,399,308,30,455,154,39,342,492,40,441,290,252,233,443,80,315,152,36,499,428,29,327,467,37,199,232,34,149,115,126,403,172,24,415,306,34,388,70,96,224,352,293,468,465,133,481,274,33,340,253,293,394,438,67,478,427,35,475,136,284,299,219,46,472,403,476,83,84,1,106,81,390,306,132,32,115,181,10,284,319,67,118,68,480,241,221,161,387,499,77,366,163,30,189,291,48,319,171,20,493,227,20,253,428,224,464,381,211,334,255,31,235,117,29,486,112,28,92,79,367,389,220,37,236,261,52,255,510,151,220,391,59,261,237,16,454,346,287,412,495,295,290,435,27,395,449,27,95,145,11,191,194,212,207,168,35,382,322,15,113,486,400,504,505,277};

const uint8_t face_right[512*3]   = {6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,5,6,6,5,6,6,6,6,6,5,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,5,6,6,5,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,5,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,5,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5,6,6,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};

const IsomerspaceForcefield::real_t X[512*3] = {2.518678,0.000000,41.417840,-2.195991,0.000000,41.417840,-20.331577,-19.101536,30.790174,-27.433760,-20.351534,23.770323,-16.295206,-35.183656,15.530528,-11.917473,-25.277078,30.790174,19.941766,-4.723959,36.238077,4.643759,15.451901,38.283696,-18.103255,8.943453,36.238077,-20.331577,19.101536,30.790174,-28.045501,-30.370053,6.483968,-35.504879,-22.134240,1.798944,-11.196136,-38.769917,11.080431,-19.647235,-37.008434,-2.915725,4.587467,-31.034528,27.456400,7.539574,-38.134884,15.530528,22.782014,-21.703405,27.456400,34.031665,-5.336291,23.770323,25.968186,10.901162,30.790174,11.930403,25.422779,30.790174,22.782014,21.703405,27.456400,-6.890091,30.545209,27.456400,4.587467,31.034528,27.456400,-27.433760,20.351534,23.770323,-39.609314,11.673088,6.483968,-27.240693,-27.528823,-16.647310,3.970850,-41.803010,-2.915725,20.725453,-34.773786,11.080431,27.309401,-27.779320,15.530528,37.930703,9.065490,15.530528,27.309401,27.779320,15.530528,35.664241,19.278047,11.080431,-11.196136,38.769917,11.080431,-28.045501,30.370053,6.483968,-35.504879,22.134240,1.798944,-19.647235,37.008434,-2.915725,-37.631644,8.966480,-16.647310,-35.009146,19.878015,-12.197212,8.906947,-7.679621,39.840852,-30.788497,-4.983962,-28.573181,-30.788497,4.983962,-28.573181,-2.908408,-36.668197,-20.895080,29.114773,-29.659149,-7.600750,36.730654,-20.607996,-2.915725,24.668173,-27.448048,-20.895080,40.047939,11.270565,-7.600750,22.852854,35.314293,-2.915725,29.114773,29.659149,-7.600750,16.133677,-33.149127,-20.895080,-23.586750,24.732625,-24.887104,-19.563004,-4.952882,-37.354859,-19.563004,4.952882,-37.354859,41.609699,-0.000000,6.483968,-8.176399,-22.855379,-34.846504,-16.384772,-11.824000,-37.354859,5.221618,-19.697069,-37.354859,17.987553,-16.556344,-34.846504,24.668173,27.448048,-20.895080,28.156562,19.797536,-24.887104,-2.908408,36.668197,-20.895080,-16.384772,11.824000,-37.354859,5.221618,19.697069,-37.354859,18.337896,9.121067,-37.354859,3.632269,6.127601,-42.006746,-41.814887,0.000000,1.798944,-19.956700,-24.041565,27.456400,-23.639762,-24.681614,23.770323,0.474901,40.398026,11.080431,-11.057353,-11.532309,38.283696,-15.111164,-13.428754,36.238077,-9.846226,-17.703940,36.238077,-3.404963,-6.072584,40.889965,7.203702,0.000000,40.889965,34.449456,-0.000000,-24.887104,11.445350,-2.851904,39.840852,16.121198,2.033805,38.283696,1.340011,2.041512,41.417840,25.968186,-10.901162,30.790174,-8.176399,22.855379,-34.846504,-1.017324,2.041512,41.417840,-12.440615,-31.888326,-24.887104,-12.681020,-20.663010,33.729722,-2.746717,11.269665,39.840852,2.788327,11.338480,39.840852,-11.057353,11.532309,38.283696,-16.135179,18.064044,33.729722,-6.881015,0.000000,40.889965,-10.984768,3.350576,39.840852,-20.175346,0.000000,36.238077,-23.686461,-14.700705,30.790174,-22.138526,-9.726363,33.729722,-27.122963,-15.436941,27.456400,-15.455655,3.868160,38.283696,-32.306695,-11.022756,23.770323,11.445350,2.851904,39.840852,-30.139897,-20.876181,19.778298,4.643759,-15.451901,38.283696,-30.286077,-26.552841,11.080431,13.762699,39.153152,-7.600750,-22.379855,-33.526007,11.080431,-20.754945,-36.393857,1.798944,-15.455655,-3.868160,38.283696,-30.703813,-5.484033,27.456400,-8.790543,-37.796447,15.530528,0.125854,-34.288094,23.770323,-3.499564,-27.774561,30.790174,3.632298,-20.038298,36.238077,11.848603,-11.057183,38.283696,-24.086830,-33.615358,6.483968,11.888506,-16.614890,36.238077,-23.586750,-24.732625,-24.887104,17.982127,-9.797989,36.238077,20.444139,-19.324510,30.790174,24.344514,-2.657174,33.729722,16.121198,-2.033805,38.283696,-34.988807,-16.527994,15.530528,27.637793,5.465637,30.790174,31.384198,2.804613,27.456400,2.788327,-11.338480,39.840852,11.848603,11.057183,38.283696,-35.009146,-19.878015,-12.197212,11.888506,16.614890,36.238077,21.035032,12.497017,33.729722,14.618888,19.566954,33.729722,-17.469715,-37.511491,6.483968,-2.067212,20.214214,36.238077,-34.126769,0.000000,23.770323,3.632298,20.038298,36.238077,40.047939,-11.270565,-7.600750,-3.499564,27.774561,30.790174,-15.292602,27.274678,27.456400,-14.156905,31.155456,23.770323,-17.549510,32.253771,19.778298,-34.886889,11.207206,19.778298,-23.639762,24.681614,23.770323,3.728808,6.071903,40.889965,-22.385874,29.079254,19.778298,-22.138526,9.726363,33.729722,-23.686461,14.700705,30.790174,-25.379810,29.263581,15.530528,-2.067212,-20.214214,36.238077,-30.703813,5.484033,27.456400,-34.886889,-11.207206,19.778298,-27.122963,15.436941,27.456400,-36.635125,0.000000,19.778298,-39.791737,-5.987502,11.080431,-38.252027,5.755055,15.530528,3.970850,41.803010,-2.915725,-9.583172,-35.482735,19.778298,-41.287012,0.000000,6.483968,-40.188857,-11.570017,1.798944,37.930703,-9.065490,15.530528,-40.940692,-5.346858,-7.600750,-25.552349,11.119104,-31.906955,-38.987732,-9.972397,-12.197212,-32.030503,-26.938615,-2.915725,-32.428933,-25.609374,-7.600750,22.852854,-35.314293,-2.915725,-10.721440,-39.994139,6.483968,-22.385874,-29.079254,19.778298,-14.957245,-39.159063,1.798944,7.539574,38.134884,15.530528,-4.780440,-41.152703,6.483968,20.444139,19.324510,30.790174,-3.057398,-41.852642,1.798944,-18.103255,-8.943453,36.238077,-33.919163,18.634021,-16.647310,-27.327120,20.495343,-24.887104,-18.840482,35.651500,-12.197212,-15.121666,23.478886,-31.906955,-27.739415,-2.525110,-31.906955,6.122141,23.587183,33.729722,-32.030503,26.938615,-2.915725,-20.051450,2.241408,-37.354859,-13.158648,-39.806806,-2.915725,-23.437890,15.096511,-31.906955,38.655173,-15.368515,6.483968,13.762699,-39.153152,-7.600750,-17.295237,-26.038439,-28.573181,-21.248022,-11.555317,-34.846504,1.340011,-2.041512,41.417840,-38.175444,-17.096628,-2.915725,-15.121666,-23.478886,-31.906955,-19.337238,36.575559,-7.600750,-5.678932,-30.799735,-28.573181,5.068118,-30.962172,-28.573181,-2.921235,-27.844680,-31.906955,24.490057,-0.000000,-34.846504,24.344514,2.657174,33.729722,18.653673,-21.044292,-31.906955,-23.684887,-4.821162,-34.846504,31.215162,-25.841037,-12.197212,8.906947,7.679621,39.840852,20.725453,34.773786,11.080431,22.160310,-22.333338,-28.573181,35.615608,-21.470009,-7.600750,17.982127,9.797989,36.238077,2.694780,-24.196446,-34.846504,15.583578,4.583446,-39.400477,-40.940692,5.346858,-7.600750,-39.609314,-11.673088,6.483968,40.144979,5.779936,-12.197212,27.689707,-5.197852,-31.906955,-12.681020,20.663010,33.729722,-10.984768,-3.350576,39.840852,14.061854,39.607825,-2.915725,-37.631644,-8.966480,-16.647310,-3.057398,41.852642,1.798944,-5.110007,15.200861,38.283696,41.609699,-0.000000,-7.600750,-38.114800,15.902924,-7.600750,-40.237899,0.000000,-12.197212,35.332473,-10.815347,-20.895080,-15.473726,3.794460,-39.400477,-7.343205,40.763316,-7.600750,35.225088,16.710520,-16.647310,29.376384,25.596661,-16.647310,6.122141,-23.587183,33.729722,-27.350782,5.283126,30.790174,8.925360,39.437176,-12.197212,-22.798768,31.329554,-16.647310,2.375212,27.927178,30.790174,-20.468454,30.469518,-20.895080,7.672093,36.021781,-20.895080,40.144979,-5.779936,-12.197212,-32.306695,11.022756,23.770323,5.822671,41.059902,6.483968,29.270075,28.013935,11.080431,35.664241,-19.278047,11.080431,-24.360155,-33.416497,-7.600750,21.816725,26.584189,-24.887104,-20.051450,-2.241408,-37.354859,-39.791737,5.987502,11.080431,7.290968,-27.092378,-31.906955,-28.441889,30.722290,1.798944,28.032936,19.971205,23.770323,36.879034,2.406504,19.778298,27.637793,-5.465637,30.790174,18.653673,21.044292,-31.906955,-29.280792,-25.335121,15.530528,-5.678932,30.799735,-28.573181,-9.583172,35.482735,19.778298,-32.295879,11.054566,-24.887104,-1.617222,40.360073,-12.197212,-28.441889,-30.722290,1.798944,-1.617222,-40.360073,-12.197212,-27.350782,-5.283126,30.790174,-9.958336,5.749283,-40.957633,35.615608,21.470009,-7.600750,-9.958336,-5.749283,-40.957633,-18.529277,-15.573920,-34.846504,10.612146,37.409737,-16.647310,-7.343205,-40.763316,-7.600750,40.760437,10.663845,-2.915725,-6.881015,0.000000,-42.006746,38.812702,3.844499,15.530528,-19.956700,24.041565,27.456400,7.290968,27.092378,-31.906955,-5.810238,-33.764107,23.770323,10.177252,5.928216,-40.957633,38.564595,5.822194,-16.647310,34.752395,-12.547478,19.778298,-15.292602,-27.274678,27.456400,22.160310,22.333338,-28.573181,-28.017979,13.735294,-28.573181,-15.812424,1.921474,-39.400477,0.313599,38.841788,-16.647310,-17.549510,-32.253771,19.778298,-34.988807,16.527994,15.530528,14.022443,14.881225,-37.354859,22.547915,9.525110,-34.846504,-16.135179,-18.064044,33.729722,-18.170952,-36.000358,11.080431,-10.211129,12.298989,-39.400477,-14.957245,39.159063,1.798944,14.022443,-14.881225,-37.354859,17.522889,-29.567743,-24.887104,0.474901,-40.398026,11.080431,-2.951765,-38.717131,15.530528,24.227927,14.339736,-31.906955,26.306685,32.839383,1.798944,0.166539,-11.638820,-40.957633,-25.379810,-29.263581,15.530528,-25.552349,-11.119104,-31.906955,-5.810238,33.764107,23.770323,5.865389,10.145246,-40.957633,-19.337238,-36.575559,-7.600750,11.993163,10.902356,-39.400477,-27.739415,2.525110,-31.906955,14.595692,-24.009958,-31.906955,7.623338,-41.307658,1.798944,29.376384,-25.596661,-16.647310,-6.934262,-41.372169,-2.915725,35.332473,10.815347,-20.895080,11.006291,17.203721,-37.354859,11.930403,-25.422779,30.790174,38.564595,-5.822194,-16.647310,-10.211129,-12.298989,-39.400477,3.632269,-6.127601,-42.006746,22.480026,31.789685,-16.647310,5.068118,30.962172,-28.573181,-4.733526,23.831210,33.729722,20.498032,0.000000,-37.354859,-18.529277,15.573920,-34.846504,0.166539,11.638820,-40.957633,-14.599564,-33.706018,-20.895080,-24.086830,33.615358,6.483968,32.746226,10.672398,-24.887104,16.180942,38.227459,6.483968,27.689707,5.197852,-31.906955,7.203702,0.000000,-42.006746,11.993163,-10.902356,-39.400477,-3.362916,6.097082,-42.006746,-32.295879,-11.054566,-24.887104,-34.054746,2.221235,-24.887104,-34.054746,-2.221235,-24.887104,0.874120,16.073122,-39.400477,18.647397,34.160994,-16.647310,-15.111164,13.428754,36.238077,12.239127,34.757836,19.778298,35.225088,-16.710520,-16.647310,14.618888,-19.566954,33.729722,25.568093,-32.748484,6.483968,-5.296170,-19.590724,-37.354859,-38.114800,-15.902924,-7.600750,27.736610,-14.910975,-28.573181,7.376435,14.380395,-39.400477,-41.406587,5.840471,-2.915725,-17.295237,26.038439,-28.573181,21.993426,29.619929,19.778298,22.547915,-9.525110,-34.846504,10.177252,-5.928216,-40.957633,30.876599,20.262110,-20.895080,31.059575,5.294508,-28.573181,28.558076,13.279987,27.456400,16.133677,33.149127,-20.895080,-23.827351,-34.446284,-2.915725,5.822671,-41.059902,6.483968,-22.151487,22.019766,-28.573181,7.376435,-14.380395,-39.400477,31.215162,25.841037,-12.197212,11.006291,-17.203721,-37.354859,-36.441799,-3.766964,-20.895080,-5.296170,19.590724,-37.354859,-28.017979,-13.735294,-28.573181,24.208383,32.462882,-12.197212,27.736610,14.910975,-28.573181,29.814219,-21.787314,19.778298,36.947472,-16.699688,-12.197212,17.700535,-36.393345,-12.197212,-11.277414,-37.119570,-16.647310,18.647397,-34.160994,-16.647310,31.059575,-5.294508,-28.573181,2.518678,0.000000,-42.534622,21.880859,35.301966,-7.600750,18.337896,-9.121067,-37.354859,24.227927,-14.339736,-31.906955,15.583578,-4.583446,-39.400477,-14.599564,33.706018,-20.895080,18.056731,34.474089,15.530528,-11.917473,25.277078,30.790174,-3.362916,-6.097082,-42.006746,-15.473726,-3.794460,-39.400477,-23.437890,-15.096511,-31.906955,41.793465,5.363803,1.798944,-33.919163,-18.634021,-16.647310,10.318140,32.749262,23.770323,-9.846226,17.703940,36.238077,32.746226,-10.672398,-24.887104,29.270075,-28.013935,11.080431,5.865389,-10.145246,-40.957633,9.744773,-32.921612,-24.887104,-14.156905,-31.155456,23.770323,-23.684887,4.821162,-34.846504,-12.440615,31.888326,-24.887104,21.035032,-12.497017,33.729722,2.694780,24.196446,-34.846504,-11.477478,0.000000,-40.957633,-2.921235,27.844680,-31.906955,-15.812424,-1.921474,-39.400477,-11.277414,37.119570,-16.647310,-34.262322,12.999667,-20.895080,14.061854,-39.607825,-2.915725,-28.306015,-28.665456,-12.197212,-24.360155,33.416497,-7.600750,-32.428933,25.609374,-7.600750,-27.240693,27.528823,-16.647310,34.128842,-18.839233,15.530528,-1.017324,-2.041512,41.417840,-23.827351,34.446284,-2.915725,-6.890091,-30.545209,27.456400,-6.934262,41.372169,-2.915725,-4.780440,41.152703,6.483968,34.128842,18.839233,15.530528,17.530485,38.214093,1.798944,14.255127,-28.001746,27.456400,10.380776,22.078258,-34.846504,10.841396,38.961972,11.080431,0.874120,-16.073122,-39.400477,30.494821,29.015239,-2.915725,34.031665,5.336291,23.770323,9.744773,32.921612,-24.887104,-10.721440,39.994139,6.483968,-29.280792,25.335121,15.530528,33.688430,25.257442,1.798944,33.292544,24.905616,6.483968,39.041413,15.822266,1.798944,10.380776,-22.078258,-34.846504,-21.248022,11.555317,-34.846504,-0.449974,34.282663,-24.887104,38.655173,15.368515,6.483968,13.465002,-28.385651,-28.573181,18.056731,-34.474089,15.530528,-24.167370,0.000000,33.729722,-41.406587,-5.840471,-2.915725,-30.139897,20.876181,19.778298,41.289211,5.144382,6.483968,-29.036075,22.393990,-20.895080,39.405203,-9.592615,11.080431,3.728808,-6.071903,40.889965,-2.951765,38.717131,15.530528,39.405203,9.592615,11.080431,40.760437,-10.663845,-2.915725,13.465002,28.385651,-28.573181,21.993426,-29.619929,19.778298,41.793465,-5.363803,1.798944,39.041413,-15.822266,1.798944,-8.790543,37.796447,15.530528,41.289211,-5.144382,6.483968,33.688430,-25.257442,1.798944,30.494821,-29.015239,-2.915725,2.375212,-27.927178,30.790174,33.292544,-24.905616,6.483968,26.306685,-32.839383,1.798944,24.208383,-32.462882,-12.197212,21.880859,-35.301966,-7.600750,10.318140,-32.749262,23.770323,36.730654,20.607996,-2.915725,3.797532,41.288549,-7.600750,17.522889,29.567743,-24.887104,7.672093,-36.021781,-20.895080,10.612146,-37.409737,-16.647310,22.480026,-31.789685,-16.647310,0.313599,-38.841788,-16.647310,8.925360,-39.437176,-12.197212,19.941766,4.723959,36.238077,-20.468454,-30.469518,-20.895080,-12.952287,-38.211667,-12.197212,-22.798768,-31.329554,-16.647310,1.953113,-36.752818,19.778298,-18.840482,-35.651500,-12.197212,-8.360234,-7.927476,39.840852,-27.327120,-20.495343,-24.887104,28.156562,-19.797536,-24.887104,14.595692,24.009958,-31.906955,-29.036075,-22.393990,-20.895080,-22.151487,-22.019766,-28.573181,25.568093,32.748484,6.483968,-34.262322,-12.999667,-20.895080,-36.441799,3.766964,-20.895080,3.797532,-41.288549,-7.600750,-0.449974,-34.282663,-24.887104,-38.680743,0.000000,-16.647310,40.560586,-0.000000,11.080431,-38.175444,17.096628,-2.915725,-40.188857,11.570017,1.798944,-38.252027,-5.755055,15.530528,-34.887212,22.126114,6.483968,-30.286077,26.552841,11.080431,-2.746717,-11.269665,39.840852,36.957811,-0.000000,-20.895080,-36.370688,17.248463,11.080431,-20.754945,36.393857,1.798944,36.947472,16.699688,-12.197212,7.623338,41.307658,1.798944,-17.469715,37.511491,6.483968,-16.295206,35.183656,15.530528,17.700535,36.393345,-12.197212,-22.379855,33.526007,11.080431,-3.404963,6.072584,40.889965,-4.733526,-23.831210,33.729722,-18.170952,36.000358,11.080431,-34.887212,-22.126114,6.483968,1.953113,36.752818,19.778298,-5.110007,-15.200861,38.283696,0.125854,34.288094,23.770323,32.716774,-10.761907,23.770323,20.196897,27.825371,23.770323,21.816725,-26.584189,-24.887104,14.255127,28.001746,27.456400,29.814219,21.787314,19.778298,42.137574,-0.000000,-2.915725,-38.987732,9.972397,-12.197212,34.752395,12.547478,19.778298,32.716774,10.761907,23.770323,28.032936,-19.971205,23.770323,-8.360234,7.927476,39.840852,-28.306015,28.665456,-12.197212,38.812702,-3.844499,15.530528,36.879034,-2.406504,19.778298,-12.952287,38.211667,-12.197212,28.558076,-13.279987,27.456400,17.987553,16.556344,-34.846504,20.196897,-27.825371,23.770323,17.530485,-38.214093,1.798944,12.239127,-34.757836,19.778298,-36.370688,-17.248463,11.080431,30.876599,-20.262110,-20.895080,-13.158648,39.806806,-2.915725,16.180942,-38.227459,6.483968,31.384198,-2.804613,27.456400,10.841396,-38.961972,11.080431};