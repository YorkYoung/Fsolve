#pragma once

#include <exception>
#include <functional>

namespace Omega
{
	using namespace std;

	class SolveFailed : public exception
	{
	public:
		SolveFailed() : 
			exception("Solve failed!")
		{

		}
	};

	template <typename TArg, typename TValue>
	TArg FZeroBisection(function<TValue(TArg)> func, TArg begin, TArg end, size_t maxStep, TArg error)
	{
		if (begin > end)
			swap(begin, end);

		TValue valueBegin = func(begin);
		for (size_t i = 0; i < maxStep; ++i)
		{
			TArg size = (end - begin) / 2;
			TArg solve = begin + size;
			TValue value = func(solve);
			if (value == 0 || size < error)
				return solve;
			if (valueBegin * value > 0)
			{
				begin = solve;
				valueBegin = value;
			}
			else
			{
				end = solve;
			}
		}
		throw SolveFailed();
	}

	template <typename TArg, typename TValue, typename TError>
	TArg FZeroNewton(function<TValue(TArg)> func, function<TValue(TArg)> diff, TArg begin, size_t maxStep, TError error)
	{
		for (size_t i = 0; i < maxStep; ++i)
		{
			TValue solve = begin - func(begin) / diff(begin);
			if (abs(solve - begin) < error)
				return solve;
			begin = solve;
		}
		throw SolveFailed();
	}

	template <typename TArg, typename TValue, typename TError>
	TArg FZeroSecant(function<TValue(TArg)> func, TArg begin0, TArg begin1,
		size_t maxStep, TError error)
	{
		TValue value0 = func(begin0);
		TValue value1 = func(begin1);
		for (size_t i = 0; i < maxStep; ++i)
		{
			TArg solve = begin1 - value1 * (begin1 - begin0) / (value1 - value0);
			if (abs(solve - begin1) < error)
				return solve;
			begin0 = begin1;
			value0 = value1;
			begin1 = solve;
			value1 = func(solve);
		}
		throw SolveFailed();
	}

	template <typename TArg, typename TValue, typename TError>
	TArg FZeroRegula(function<TValue(TArg)> func, TArg begin0, TArg begin1,
		size_t maxStep, TError error)
	{
		TValue value0 = func(begin0);
		TValue value1 = func(begin1);
		for (size_t i = 0; i < maxStep; ++i)
		{
			TArg solve = begin1 - value1 * (begin1 - begin0) / (value1 - value0);
			if (abs(solve - begin1) < error)
				return solve;
			TValue value = func(solve);
			if (value * value1 < 0)
			{
				begin0 = begin1;
				value0 = value1;
			}
			begin1 = solve;
			value1 = func(solve);
		}
		throw SolveFailed();
	}

	template <typename T, typename TError>
	T FixedPointIteration(function<T(T)> func, T begin, size_t maxStep, TError error)
	{
		for (size_t i = 0; i < maxStep; ++i)
		{
			T value = func(begin);
			if (abs(value - begin) < error)
				return value;
			begin = value;
		}
		throw SolveFailed();
	}

	template <typename T, typename TError>
	T FixedPointSteffensen(function<T(T)> func, T begin, size_t maxStep, TError error)
	{
		T value1, value2;
		for (int i = 0; i < maxStep; ++i)
		{
			value1 = func(begin);
			value2 = func(value1);
			T value = begin - (value1 - begin) * (value1 - begin) / (value2 - 2 * value1 + begin);
			if (abs(value - begin) < error)
				return value;
			begin = value;
		}
		throw SolveFailed();
	}
}
