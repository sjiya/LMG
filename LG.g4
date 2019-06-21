grammar LG;

WS: [ \t\r\n]+ -> skip;

ADD: '+';
SUB: '-';
MUL: '*';
DIV: '/';
AMPERSENT :'&';
BACKSLASHES : '\\\\';
EOL: '\r'?'\n';

L_PAREN: '(';
R_PAREN: ')';
L_BRACE: '{';
R_BRACE: '}';
L_BRACKET: '[';
R_BRACKET: ']';
LR_BRACE: '{}';
BAR: '|';

BEGINMATRIX: '\\begin{bmatrix}' | '\\begin{pmatrix}' | '\\begin{vmatrix}';
ENDMATRIX: '\\end{bmatrix}' | '\\end{pmatrix}' | '\\end{vmatrix}';

BEGINCASES: '\\begin{cases}';
ENDCASES: '\\end{cases}';

FUNC_LIM:  '\\lim';
LIM_APPROACH_SYM: '\\to' | '\\rightarrow' | '\\Rightarrow' | '\\longrightarrow' | '\\Longrightarrow';
FUNC_INT:  '\\int';
FUNC_SUM:  '\\sum';
FUNC_PROD: '\\prod';

FUNC_OVERLINE:'\\overline';
FUNC_XRIGHTARROW:'\\xrightarrow' | '\\xRightarrow';
FUNC_UNDERBRACE: '\\underbrace';

FUNC_LOG:  '\\log';
FUNC_LN:   '\\ln';
FUNC_SIN:  '\\sin';
FUNC_COS:  '\\cos';
FUNC_TAN:  '\\tan';
FUNC_CSC:  '\\csc';
FUNC_SEC:  '\\sec';
FUNC_COT:  '\\cot';

FUNC_ARCSIN: '\\arcsin';
FUNC_ARCCOS: '\\arccos';
FUNC_ARCTAN: '\\arctan';
FUNC_ARCCSC: '\\arccsc';
FUNC_ARCSEC: '\\arcsec';
FUNC_ARCCOT: '\\arccot';

FUNC_SINH: '\\sinh';
FUNC_COSH: '\\cosh';
FUNC_TANH: '\\tanh';
FUNC_ARSINH: '\\arsinh';
FUNC_ARCOSH: '\\arcosh';
FUNC_ARTANH: '\\artanh';


FUNC_SQRT: '\\sqrt';


FUNC_LFLOOR:'\\lfloor';
FUNC_RFLOOR:'\\rfloor';
FUNC_BIG: '\\Big';

CMD_TIMES: '\\times';
CMD_CDOT:  '\\cdot';
CMD_DIV:   '\\div';
CMD_FRAC:  '\\frac';
CMD_CFRAC: '\\cfrac';

CMD_MATHIT: '\\mathit' | '\\mathrm' | '\\mathbf' | '\\mathfrak' | '\\mathsf' | '\\mathtt' | '\\mathcal' | '\\text' |'\\mbox' |'\\textit' | '\\textbf';
STAR: '\u002A';
UNDERSCORE: '_';
CARET: '^';
COLON: ':';

fragment WS_CHAR: [ \t\r\n];
DIFFERENTIAL: 'd' WS_CHAR*? ([a-zA-Z] | '\\' [a-zA-Z]+);

LETTER: [a-zA-Z];
fragment DIGIT: [0-9];
NUMBER:
    DIGIT+ (',' DIGIT DIGIT DIGIT)*
    | DIGIT* (',' DIGIT DIGIT DIGIT)* '.' DIGIT+;
//TEXT: ('a..z'|'A..Z'|'0..9'|','| '\n' |' ')+;
//WORD:   [A-Za-z0-9' ]+;
INVTIMES: [a-zA-Z0-9]+?;
EQUAL: '=';
LT: '<';
LTE: '\\leq';
GT: '>';
GTE: '\\geq';
NE: '\\neq';
IN: '\\in';

BANG: '!';

SYMBOL: '\\' [a-zA-Z]+;

math: relation;

relation:
    relation (EQUAL | LT | LTE | GT | GTE |NE |IN |LIM_APPROACH_SYM) relation
    | expr;

equality:
    expr EQUAL expr;

expr: additive;

additive:
    additive (ADD | SUB) additive
    | mp;

// mult part
mp:
    mp (MUL | CMD_TIMES | CMD_CDOT | DIV | CMD_DIV | COLON) mp
    | unary;

mp_nofunc:
    mp_nofunc (MUL | CMD_TIMES | CMD_CDOT | DIV | CMD_DIV | COLON) mp_nofunc
    | unary_nofunc;

unary:
    (ADD | SUB) unary
    | postfix+;

unary_nofunc:
    (ADD | SUB) unary_nofunc
    | postfix postfix_nofunc*;

postfix: exp postfix_op*;
postfix_nofunc: exp_nofunc postfix_op*;
postfix_op: BANG | eval_at;

eval_at:
    BAR (eval_at_sup | eval_at_sub | eval_at_sup eval_at_sub);

eval_at_sub:
    UNDERSCORE L_BRACE
    (expr | equality)
    R_BRACE;

eval_at_sup:
    CARET L_BRACE
    (expr | equality)
    R_BRACE;

exp:
    exp CARET (atom | L_BRACE expr R_BRACE) subexpr?
    | comp;

exp_nofunc:
    exp_nofunc CARET (atom | L_BRACE expr R_BRACE) subexpr?
    | comp_nofunc;

comp:
    group
    | abs_group
    | func
    | atom
    | frac
    | cfrac;

comp_nofunc:
    group
    | abs_group
    | atom
    | frac
    | cfrac;

group:
    L_PAREN expr R_PAREN 
    | L_BRACKET expr R_BRACKET
    | L_BRACE expr R_BRACE;

abs_group: BAR expr BAR;

atom: (LETTER | SYMBOL) subexpr? | NUMBER | DIFFERENTIAL | mathit |','|'*'|INVTIMES; //| (LETTER | SYMBOL) supexpr?

mathit: CMD_MATHIT L_BRACE mathit_text R_BRACE;
mathit_text: LETTER*;

frac:
    CMD_FRAC L_BRACE
    upper=expr*
    R_BRACE L_BRACE
    lower=expr*
    R_BRACE;
cfrac:
    CMD_CFRAC L_BRACE
    upper=expr*
    R_BRACE L_BRACE
    lower=expr*
    R_BRACE;

matrix : matrixdef;
matrixdef : BEGINMATRIX matrixrows ENDMATRIX (BEGINMATRIX matrixrows ENDMATRIX)*;
matrixrows : matrixrow (rowbreak matrixrow)*;
matrixrow : relation (columnbreak relation )*;

casesdef: BEGINCASES caserows ENDCASES;
caserows: caserow (rowbreak caserow)*;
caserow: relation (columnbreak relation)*;

rowbreak :  (BACKSLASHES|EOL);
columnbreak : AMPERSENT;

chemele: LR_BRACE supexpr subexpr mathit;

func_normal:
    FUNC_LOG | FUNC_LN
    | FUNC_SIN | FUNC_COS | FUNC_TAN
    | FUNC_CSC | FUNC_SEC | FUNC_COT
    | FUNC_ARCSIN | FUNC_ARCCOS | FUNC_ARCTAN
    | FUNC_ARCCSC | FUNC_ARCSEC | FUNC_ARCCOT
    | FUNC_SINH | FUNC_COSH | FUNC_TANH
    | FUNC_ARSINH | FUNC_ARCOSH | FUNC_ARTANH;

func_cmd:
	FUNC_LFLOOR | FUNC_RFLOOR;
func_bracket:
	FUNC_BIG;

func:
    func_normal
    (subexpr? supexpr? | supexpr? subexpr?)
    (L_PAREN func_arg R_PAREN | func_arg_noparens)

    | (LETTER | SYMBOL) subexpr? // e.g. f(x)
    L_PAREN args R_PAREN

    | FUNC_INT
    (subexpr supexpr | supexpr subexpr)?
    (additive? DIFFERENTIAL | frac |cfrac | additive)

    | FUNC_SQRT
    (L_BRACKET root=expr R_BRACKET)?
    L_BRACE base=expr R_BRACE
	
    | (FUNC_SUM | FUNC_PROD)
    (subeq supexpr | supexpr subeq)
    mp
    | FUNC_LIM limit_sub mp

    | matrix
    | FUNC_OVERLINE L_BRACE expr R_BRACE
    | casesdef
    | FUNC_XRIGHTARROW L_BRACE expr R_BRACE
    | chemele
    | FUNC_UNDERBRACE L_BRACE relation (relation)* R_BRACE subexpr
    
    |func_bracket func_cmd relation (relation)* func_bracket func_cmd
    | func_cmd 
    | func_bracket | func_cmd |(L_PAREN expr | L_BRACE | L_BRACKET) additive? | func_cmd | func_bracket (R_PAREN | R_BRACE | R_BRACKET);
args: (expr ',' args) | expr;

limit_sub:
    UNDERSCORE L_BRACE
    (LETTER | SYMBOL)
    LIM_APPROACH_SYM
    expr (CARET L_BRACE (ADD | SUB) R_BRACE)?
    R_BRACE;

func_arg: expr | (expr ',' func_arg);
func_arg_noparens: mp_nofunc;

subexpr: UNDERSCORE (atom | L_BRACE expr R_BRACE);
supexpr: CARET (atom | L_BRACE expr R_BRACE);

subeq: UNDERSCORE L_BRACE equality R_BRACE;
supeq: UNDERSCORE L_BRACE equality R_BRACE;
