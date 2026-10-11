//
// Program.h from Evaluator
// 

#pragma once

#include <cstdint>
#include <vector>
#include <stack>
#include <random>
#include "vessl/vessl.h"

namespace Program
{
	// type of the string expression for Compile
	typedef char	 Char;
	// type of the value returned by evaluation
	typedef uint64_t Value;

	enum CompileError
	{
		CE_NONE,
		CE_MISSING_PAREN,
		CE_MISSING_BRACKET,
		CE_MISSING_BRACE,
		CE_MISSING_COLON_IN_TERNARY,
		CE_UNEXPECTED_CHAR,
		CE_FAILED_TO_PARSE_NUMBER,
		CE_ILLEGAL_ASSIGNMENT,
		CE_ILLEGAL_STATEMENT_TERMINATION, // found a semi-colon where one isn't allowed
		CE_ILLEGAL_VARIABLE_NAME, // found an uppercase letter where we expected a lowercase one
		CE_MISSING_PUT, // the program does not contain any PUT instructions
		CE_TOO_MANY_OPERATIONS, // the compiled program contained more operations that maximum amount specific by the user
	};

	static const char * GetErrorString(CompileError error);

	struct RuntimeError
	{
		enum Code
		{
			RE_NONE,
			RE_DIVIDE_BY_ZERO,
			RE_MISSING_OPERAND, // the stack did not contain data for the operation
			RE_MISSING_OPCODE, // the implementation for an Op::Code is missing
			RE_INCONSISTENT_STACK, // the stack was not empty after popping the final result
			RE_EMPTY_PROGRAM, // the program has no instructions to execute
			RE_GET_OUT_OF_BOUNDS, // index for GET was bigger than the size of the results array
			RE_PUT_OUT_OF_BOUNDS, // index for PUT was bigger than the size of the results array
		};

		Code code = RE_NONE;
		Value value;

		operator bool() {return code != RE_NONE; }
	};

	static const char * GetErrorString(RuntimeError error);

	struct Op
	{
	public:

		enum Code
		{
			NOP, // no operation (doesn't pop from the stack or push to it)
			PSH, // push a constant value onto the stack (eg when a numeric value is used)
			PEK, // get the value at a memory address
			POK, // set the value at a memory address
			FRQ,
			SQR,
			SIN,
			TRI,
			NEG,
			MUL,
			DIV,
			MOD,
			ADD,
			SUB,
			BSL,
			BSR,
			AND,
			OR,
			XOR,
			CEQ, // compare ==
			CNE, // compare !=
			CLT, // compare <
			CLE, // compare <=
			CGT, // compare >
			CGE, // compare >=
			CND, // conditional - ?: and ?
			POP, // ; (pop a value from the stack and do nothing with it)
			GET, // get the the current value of a result. eg [0] or [1].
			PUT, // assign to an output result using [0] = expression.
			RND, // random number operator - operand is used to wrap value returned by rand() - so like Random.Range(0, operand)
			CCV, // use the operand to look up the current value of a midi control change value, eg 'a = C1'
			VCV, // use the operand to look up the current value of a "voltage" control value, eg 'a = V5'
			NOT,
			COM,
			JMP, // JMP to the address indicated by val
		};

		// need default constructor or we can't use vector
		Op() : code(PSH), val(0) {}
		Op(Code _code, Value _val) : code(_code), val(_val) {}

		Code code;
		// this is mutable because the compiler needs to be able to change it in some cases
		Value val;
	};
	
	// fixed-size stack
	template<typename T, size_t N>
	class Stack
	{
	public:
		VESSL_INLINE bool empty() const { return t_==0; }
		VESSL_INLINE size_t size() const { return t_; }
		VESSL_INLINE void push(const T& v) { s_[t_++] = v; }
		VESSL_INLINE void pop() { VASSERT(!empty(), "pop on an empty stack"); --t_; }
		VESSL_INLINE const T& top() const { VASSERT(!empty(), "top on an empty stack"); return s_[t_-1]; }
		VESSL_INLINE void clear() { t_ = 0; }

	private:
		T s_[N];
		size_t t_ = 0;
	};

	// returns number of operations in compiled program
	template<size_t OpsMax, size_t MemSize>
	static size_t Compile(const Char* source, std::vector<Op>& outOps, CompileError& outError, int& outErrorPosition);

	// get the address in memory of a variable declared in a program with a particular userMemorySize.
	static Value GetAddress(const Char var, size_t userMemorySize);

	// get human-readable descriptions of errors

	template<size_t OpsMax, size_t MemSize>
	class Executable
	{
	public:
		Executable();

		bool Load(const std::vector<Op>& inOps, const size_t seed);

		uint64_t GetInstructionCount() const { return ops_size; }

		// run the program placing the value it evaluates to into the results array.
		// count is provided so that we can prevent the program from overrunning the array.
		RuntimeError Run(Value* results, const size_t size);

		// get the current value of a var, eg Get('t')
		Value Get(const Char var) const;
		// set the value of a var, eg Set('m', 128)
		void  Set(const Char var, const Value value);

		Value& operator[](const Char var);

		// get the value at this memory address
		Value Peek(const Value address) const;
		// set the value at this memory address
		void  Poke(const Value address, const Value value);

		// get/set a control change. these are accessible in the program with the 'C' operator
		Value GetCC(const Value idx) const;
		void  SetCC(const Value idx, const Value value);

		// get/set a "voltage" change. these are accessible in the program with the 'V' operator
		Value GetVC(const Value idx) const;
		void  SetVC(const Value idx, const Value value);

	private:
		RuntimeError Exec(const Op& op, Value* results, size_t size);

		static constexpr size_t kCCSize = 128;
		static constexpr size_t kVCSize = 8;
		static constexpr size_t kStackMax = 32;

		static constexpr size_t userMemSize = MemSize; // how much of mem is "user" memory
		static constexpr size_t memSize = MemSize + 256; // the actual size of mem

		// the memory space - read/write memory for the program (use Peek/Poke from C++)
		// this includes "user" memory accessible with @, where @0 maps to mem[0]
		// and also includes "variable" memory accessible with lowercase letters like 'a', 'b', 'c', etc.
		// it is also possible to access variable values with @ if you know the address of the variable.
		// for safety, we always wrap the address to the size of the array to prevent invalid access.
		Value mem[MemSize + 256];
		// memory for storing MIDI CC values - readonly from within a program
		Value cc[kCCSize];
		// memory for storing VC values = readonly from within a program
		Value vc[kVCSize];
		// the compiled code
		Op ops[OpsMax];
		// the execution stack (reused each time Run is called)
		Stack<Value, kStackMax> stack;

		size_t ops_size; // how many ops we actuall have
		size_t pc; // program counter, stored here because it can be changed by TRN and JMP
		float  fc; // pre-computed coefficent for FRQ

		// rng because rand() doesn't generate a large enough range
		std::default_random_engine rng;
	};
};

#include "Program.inl"