module list_mod
	use kinds_mod
	implicit none

	type::list_element
		type(list_element),pointer::left => null()
		type(list_element),pointer::right => null()						
		real(wp)::value = 0.0_wp
	end type

	type::list
		integer::length
		type(list_element),private,pointer::first => null()
		type(list_element),private,pointer::last => null()
		type(list_element),private,pointer::ptr => null()
		integer::index = 0
	contains
		procedure::new => list_construct
		procedure::fill => list_fill
		procedure::finalize => list_finalize
		procedure::traverse => list_traverse
		procedure::insert => list_insert
		procedure::remove => list_remove
		procedure::get => list_get
		procedure::set => list_set
		procedure::push => list_push
		procedure::freeze => list_freeze
	end type

contains

	subroutine list_freeze(this,frozen)
		class(list),intent(inout)::this
		real(wp),allocatable,dimension(:)::frozen
		integer::i
		allocate(frozen(this%length))
		
		do i=1,this%length
			frozen(i) = this%get(i)
		end do
	end subroutine

	function list_traverse(this,index) result(v)
		class(list),intent(inout)::this
		integer,intent(in)::index
		integer::d
		type(list_element),pointer::v

		if(index > this%length .or. index < 1) then
			error stop 'Index out of bounds.'
		end if

		d = abs(this%index-index)

		if(d > index - 1 .and. index < this%length/2) then
			this%ptr => this%first
			this%index = 1
		else if (d > this%length - index) then
			this%ptr => this%last
			this%index = this%length
		end if

		call list_traverse_to(this,this%index,this%ptr,index)
		v => this%ptr
	end function
	
	subroutine list_traverse_to(this,index,e,to)
		class(list),intent(in)::this
		integer,intent(inout)::index
		type(list_element),pointer,intent(inout)::e
		integer,intent(in)::to

		if(index <= to) then
			do while(.true.)
				if(index >= to) then
					exit
				end if
				e => e%right
				index = index + 1
			end do
		else
			do while(.true.)
				if(index <= to) then
					exit
				end if
				e => e%left
				index = index - 1
			end do		
		endif
	end subroutine

	function list_get(this, index) result(v)
		class(list)::this
		integer, intent(in)::index
		real(wp)::v
!		v => this%traverse(index)%value
		type(list_element),pointer::e
		e => this%traverse(index)
		v = e%value
	end function
	
	subroutine list_set(this, index, value)
		class(list)::this
		integer, intent(in)::index
		real(wp), intent(in)::value
		type(list_element),pointer::e
		e => this%traverse(index)
		e%value = value
	end subroutine

	subroutine list_insert(this, index, value) ! inserts after index
		class(list)::this
		integer, intent(in)::index
		real(wp), intent(in)::value
		type(list_element), pointer::left, inserted, right

		allocate(inserted)

		if(index==0) then
			this%first%left=>inserted
			inserted%right => this%first
			this%first => inserted
		else if(index==this%length) then
			this%last%right => inserted
			inserted%left => this%last
			this%last => inserted
		else
			left => this%traverse(index)
			right => left%right
			left%right => inserted
			right%left => inserted
			inserted%left => left
			inserted%right => right
		end if
		inserted%value = value
		
		this%length = this%length + 1
	end subroutine

	subroutine list_remove(this, index)
		class(list)::this
		integer, intent(in)::index
		type(list_element), pointer::left, removed, right

		if(this%length <= 2) then
			error stop 'Minimum number of elements in list is 2'
		end if

		removed => this%traverse(index)
		
		if(index==1) then
			this%first => removed%right
			nullify(this%first%left)
		else if(index==this%length) then
			this%last => removed%left
			nullify(this%last%right)
		else
			left => removed%left
			right => removed%right

			left%right => right
			right%left => left
		end if

		deallocate(removed)
		this%length = this%length - 1
	end subroutine

	subroutine list_finalize(this)
		class(list)::this
		integer::i
		type(list_element), pointer:: deleted, right

		right => this%first
		do i = 1, this%length
			deleted => right
			if(i/=this%length) then
				right => deleted%right
			end if
			deallocate(deleted)
		end do
	end subroutine

	subroutine list_construct(this)
		class(list)::this
		this%length = 2
		allocate(this%first)
		allocate(this%last)

		this%first%right => this%last
		this%last%left => this%first
		this%ptr => this%first
		this%index = 1
	end subroutine
	
	subroutine list_push(this,value)
		class(list)::this
		real(wp),intent(in)::value
		
		call this%insert(this%length,value)
	end subroutine
	
	subroutine list_fill(this,array)
		class(list)::this
		real(wp),dimension(:),intent(in)::array
		integer::i
		
		if(associated(this%first)) then
			call this%finalize()
		end if
		call this%new()
		
		call this%set(1,array(1))
		call this%set(2,array(2))
		do i = 3, size(array)
			call this%push(array(i))
		end do
	end subroutine
end module
