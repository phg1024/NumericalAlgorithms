function lessThan(a)
	comp(b) = b<a?true:false
	return comp
end

function greaterThan(a)
	comp(b) = b>a?true:false
	return comp
end

function equal(a)
	comp(b) = a==b
	return comp
end

function merge(A, B)
	if length(A) == 0
		return B
	elseif length(B) == 0
		return A
	else
		return [A B]
	end
end

function merge(A, B, C)
	return merge(merge(A, B), C)
end

function quicksort(A)
	if length(A) <= 1
		return A
	end

	pivot = A[1]
	left = filter(lessThan(pivot), A)
	right = filter(greaterThan(pivot), A)
	middle = filter(equal(pivot), A)
	return merge(quicksort(left), middle, quicksort(right))
end
